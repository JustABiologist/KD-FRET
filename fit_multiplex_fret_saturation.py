#!/usr/bin/env python
"""Fit FRET saturation curves for all worksheets in a multiplex results workbook.

The fitted model is:

    Y = Bmax * X / (Kd + X)

Here Kd is in the same units as X, which is AccBB_BG for both fitted plots.

Two analyses are produced:
1. Raw: x = AccBB_BG, y = CorrFRET.
2. Normalized response: x = AccBB_BG, y = CorrFRET / DonBB_BG.
"""

from __future__ import annotations

import argparse
import ctypes
import math
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


DEFAULT_WORKBOOK = Path(
    r"C:\Users\gruenf\Desktop\KD-FRET\PROCESSED\04052026_Test_poolfix\multiplex_results.xlsx"
)
DEFAULT_OUTPUT_DIR = Path("analysis_outputs")
REQUIRED_COLUMNS = (
    "analysis_included",
    "in_bleach_roi",
    "AccBB_BG",
    "CorrFRET",
    "DonBB_BG",
)
BARPLOT_COLUMNS = ("analysis_included", "well", "DonBB_BG")
COLORS = {
    "Hepta-Lisa": "#2563eb",
    "Hepta-Multi": "#f97316",
    "Penta-Lisa": "#16a34a",
    "Penta-Multi": "#7c3aed",
}
LISA_IPTG_BY_WELL = {
    "Well01": 3,
    "Well02": 10,
    "Well03": 25,
    "Well04": 50,
}


@dataclass
class FitResult:
    analysis: str
    worksheet: str
    n: int
    bmax: float
    kd: float
    ns: float
    background: float
    bmax_se: float
    kd_se: float
    ns_se: float
    background_se: float
    r_squared: float
    adjusted_r_squared: float
    rmse: float
    rss: float
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    status: str


def saturation_model(
    x: np.ndarray | float,
    bmax: float,
    kd: float,
) -> np.ndarray | float:
    return bmax * x / (kd + x)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit Bmax/Kd saturation curves for multiplex KD-FRET worksheets."
    )
    parser.add_argument(
        "--workbook",
        type=Path,
        default=DEFAULT_WORKBOOK,
        help=f"Input multiplex_results.xlsx. Default: {DEFAULT_WORKBOOK}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output folder. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--no-snapshot",
        action="store_true",
        help="Read the workbook directly instead of copying a snapshot first.",
    )
    parser.add_argument(
        "--min-x",
        type=float,
        default=0.0,
        help="Minimum x value retained for fitting and plotting. Default: 0.",
    )
    parser.add_argument(
        "--point-alpha",
        type=float,
        default=0.18,
        help="Scatter point alpha. Default: 0.18.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=9.0,
        help="Scatter point size. Default: 9.",
    )
    parser.add_argument(
        "--maxfev",
        type=int,
        default=200000,
        help="Maximum function evaluations for scipy.optimize.curve_fit. Default: 200000.",
    )
    return parser.parse_args()


def copy_workbook_allow_open(source: Path, destination: Path) -> Path:
    """Copy a workbook even when Excel keeps it open on Windows."""
    destination.parent.mkdir(parents=True, exist_ok=True)

    try:
        shutil.copyfile(source, destination)
        return destination
    except PermissionError:
        pass

    if not hasattr(ctypes, "windll"):
        raise

    generic_read = 0x80000000
    share_read = 0x00000001
    share_write = 0x00000002
    share_delete = 0x00000004
    open_existing = 3
    invalid_handle_value = ctypes.c_void_p(-1).value

    create_file = ctypes.windll.kernel32.CreateFileW
    create_file.argtypes = [
        ctypes.c_wchar_p,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_uint32,
        ctypes.c_void_p,
    ]
    create_file.restype = ctypes.c_void_p

    read_file = ctypes.windll.kernel32.ReadFile
    read_file.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_uint32),
        ctypes.c_void_p,
    ]
    read_file.restype = ctypes.c_int

    close_handle = ctypes.windll.kernel32.CloseHandle
    close_handle.argtypes = [ctypes.c_void_p]
    close_handle.restype = ctypes.c_int

    handle = create_file(
        str(source),
        generic_read,
        share_read | share_write | share_delete,
        None,
        open_existing,
        0,
        None,
    )
    if handle == invalid_handle_value:
        raise PermissionError(f"Could not open workbook, even with shared read: {source}")

    try:
        with destination.open("wb") as out:
            buffer_size = 1024 * 1024
            buffer = ctypes.create_string_buffer(buffer_size)
            bytes_read = ctypes.c_uint32(0)

            while True:
                ok = read_file(handle, buffer, buffer_size, ctypes.byref(bytes_read), None)
                if not ok:
                    raise OSError(f"Windows ReadFile failed while copying {source}")
                if bytes_read.value == 0:
                    break
                out.write(buffer.raw[: bytes_read.value])
    finally:
        close_handle(handle)

    return destination


def prepare_workbook_path(workbook: Path, output_dir: Path, no_snapshot: bool) -> Path:
    workbook = workbook.expanduser().resolve()
    if no_snapshot:
        return workbook

    snapshot = output_dir / f"{workbook.stem}_snapshot.xlsx"
    return copy_workbook_allow_open(workbook, snapshot)


def validate_columns(sheet_name: str, df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Worksheet {sheet_name!r} is missing columns: {missing}")


def validate_barplot_columns(sheet_name: str, df: pd.DataFrame) -> None:
    missing = [column for column in BARPLOT_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Worksheet {sheet_name!r} is missing columns: {missing}")


def finite_filtered_frame(df: pd.DataFrame, min_x: float) -> pd.DataFrame:
    for column in REQUIRED_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    keep = df["analysis_included"].eq(1)
    keep &= df["in_bleach_roi"].eq(1)
    keep &= np.isfinite(df["AccBB_BG"])
    keep &= np.isfinite(df["CorrFRET"])
    keep &= np.isfinite(df["DonBB_BG"])
    keep &= df["DonBB_BG"].gt(0)
    keep &= df["AccBB_BG"].ge(min_x)
    return df.loc[keep, list(REQUIRED_COLUMNS)].copy()


def well_sort_key(well: object) -> tuple[str, int, str]:
    text = str(well)
    digits = "".join(char for char in text if char.isdigit())
    prefix = "".join(char for char in text if not char.isdigit())
    number = int(digits) if digits else 10**9
    return prefix, number, text


def summarize_gfp_by_well(sheet_name: str, df: pd.DataFrame) -> pd.DataFrame:
    validate_barplot_columns(sheet_name, df)
    df = df.copy()
    df["analysis_included"] = pd.to_numeric(df["analysis_included"], errors="coerce")
    df["DonBB_BG"] = pd.to_numeric(df["DonBB_BG"], errors="coerce")
    keep = df["analysis_included"].eq(1) & np.isfinite(df["DonBB_BG"])
    summary = (
        df.loc[keep]
        .groupby("well", sort=False)["DonBB_BG"]
        .agg(["count", "mean", "std"])
        .reset_index()
    )
    summary["sem"] = summary["std"] / np.sqrt(summary["count"])
    summary["worksheet"] = sheet_name
    summary["iptg_uM"] = np.nan
    if sheet_name.endswith("-Lisa"):
        summary["iptg_uM"] = summary["well"].map(LISA_IPTG_BY_WELL)
    summary["sort_key"] = summary["well"].map(well_sort_key)
    summary = summary.sort_values("sort_key").drop(columns=["sort_key"])
    return summary[["worksheet", "well", "iptg_uM", "count", "mean", "std", "sem"]]


def initial_guess(x: np.ndarray, y: np.ndarray) -> list[float]:
    y_high = float(np.nanpercentile(y, 95))
    bmax = max(y_high, np.finfo(float).eps)

    positive_x = x[np.isfinite(x) & (x > 0)]
    kd = float(np.nanmedian(positive_x)) if positive_x.size else 1.0
    kd = max(kd, np.finfo(float).eps)

    return [bmax, kd]


def fit_curve(
    analysis: str,
    worksheet: str,
    x: np.ndarray,
    y: np.ndarray,
    maxfev: int,
) -> FitResult:
    n = int(x.size)
    if n < 5:
        nan = float("nan")
        return FitResult(
            analysis,
            worksheet,
            n,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            nan,
            float(np.nanmin(x)) if n else nan,
            float(np.nanmax(x)) if n else nan,
            float(np.nanmin(y)) if n else nan,
            float(np.nanmax(y)) if n else nan,
            "too_few_points",
        )

    try:
        popt, pcov = curve_fit(
            saturation_model,
            x,
            y,
            p0=initial_guess(x, y),
            bounds=(
                [0.0, np.finfo(float).eps],
                [np.inf, np.inf],
            ),
            maxfev=maxfev,
        )
        fitted = saturation_model(x, *popt)
        residuals = y - fitted
        rss = float(np.sum(residuals**2))
        tss = float(np.sum((y - np.mean(y)) ** 2))
        r_squared = 1.0 - (rss / tss) if tss > 0 else float("nan")
        p = len(popt)
        adjusted_r_squared = (
            1.0 - ((1.0 - r_squared) * (n - 1) / (n - p - 1))
            if n > p + 1 and np.isfinite(r_squared)
            else float("nan")
        )
        rmse = float(math.sqrt(rss / n))
        parameter_variances = np.diag(pcov)
        parameter_variances = np.where(parameter_variances >= 0, parameter_variances, np.nan)
        parameter_se = np.sqrt(parameter_variances)
        status = "ok_curve_fit_hyperbola"
    except Exception as exc:
        popt = np.full(2, np.nan)
        parameter_se = np.full(2, np.nan)
        rss = float("nan")
        r_squared = float("nan")
        adjusted_r_squared = float("nan")
        rmse = float("nan")
        status = f"fit_failed: {type(exc).__name__}: {exc}"

    return FitResult(
        analysis=analysis,
        worksheet=worksheet,
        n=n,
        bmax=float(popt[0]),
        kd=float(popt[1]),
        ns=float("nan"),
        background=float("nan"),
        bmax_se=float(parameter_se[0]),
        kd_se=float(parameter_se[1]),
        ns_se=float("nan"),
        background_se=float("nan"),
        r_squared=r_squared,
        adjusted_r_squared=adjusted_r_squared,
        rmse=rmse,
        rss=rss,
        x_min=float(np.min(x)),
        x_max=float(np.max(x)),
        y_min=float(np.min(y)),
        y_max=float(np.max(y)),
        status=status,
    )


def result_to_dict(result: FitResult) -> dict[str, object]:
    return {
        "analysis": result.analysis,
        "worksheet": result.worksheet,
        "n": result.n,
        "Bmax": result.bmax,
        "Bmax_SE": result.bmax_se,
        "Kd": result.kd,
        "Kd_SE": result.kd_se,
        "NS": result.ns,
        "NS_SE": result.ns_se,
        "Background": result.background,
        "Background_SE": result.background_se,
        "R_squared": result.r_squared,
        "Adjusted_R_squared": result.adjusted_r_squared,
        "RMSE": result.rmse,
        "RSS": result.rss,
        "x_min": result.x_min,
        "x_max": result.x_max,
        "y_min": result.y_min,
        "y_max": result.y_max,
        "status": result.status,
    }


def nice_axis_limits(values: Iterable[np.ndarray]) -> tuple[float, float]:
    combined = np.concatenate([array[np.isfinite(array)] for array in values])
    if combined.size == 0:
        return 0.0, 1.0

    low = float(np.nanpercentile(combined, 0.5))
    high = float(np.nanpercentile(combined, 99.5))
    if not np.isfinite(low) or not np.isfinite(high) or low == high:
        return float(np.nanmin(combined)), float(np.nanmax(combined))

    span = high - low
    return low - (0.05 * span), high + (0.08 * span)


def plot_analysis(
    *,
    analysis: str,
    datasets: dict[str, tuple[np.ndarray, np.ndarray]],
    results: dict[str, FitResult],
    output_path: Path,
    x_label: str,
    y_label: str,
    title: str,
    point_alpha: float,
    point_size: float,
) -> None:
    fig, ax = plt.subplots(figsize=(11.5, 8.0), constrained_layout=True)

    all_x = []
    all_y = []
    for worksheet, (x, y) in datasets.items():
        color = COLORS.get(worksheet, None)
        all_x.append(x)
        all_y.append(y)
        ax.scatter(
            x,
            y,
            s=point_size,
            alpha=point_alpha,
            color=color,
            edgecolors="none",
            rasterized=True,
            label=f"{worksheet} data (n={x.size})",
        )

        result = results[worksheet]
        if result.status.startswith("ok"):
            x_curve = np.linspace(max(0.0, result.x_min), result.x_max, 500)
            y_curve = saturation_model(
                x_curve,
                result.bmax,
                result.kd,
            )
            ax.plot(
                x_curve,
                y_curve,
                color=color,
                linewidth=2.5,
                label=f"{worksheet} fit: Kd={result.kd:.3g}, Bmax={result.bmax:.3g}",
            )

    x_low, x_high = nice_axis_limits(all_x)
    y_low, y_high = nice_axis_limits(all_y)
    ax.set_xlim(left=max(0.0, x_low), right=x_high)
    ax.set_ylim(bottom=y_low, top=y_high)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True, color="#d0d7de", linewidth=0.8, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="best", frameon=True, fontsize=9)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    svg_path = output_path.with_suffix(".svg")
    fig.savefig(svg_path)
    plt.close(fig)
    print(f"Wrote {output_path}")
    print(f"Wrote {svg_path}")


def plot_gfp_by_well(
    well_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    worksheets = list(dict.fromkeys(well_summary["worksheet"]))
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), constrained_layout=True, sharey=True)
    axes_flat = axes.ravel()

    y_max = float((well_summary["mean"] + well_summary["sem"].fillna(0)).max())
    y_limit = y_max * 1.12 if np.isfinite(y_max) and y_max > 0 else 1.0

    for axis_index, ax in enumerate(axes_flat):
        if axis_index >= len(worksheets):
            ax.set_visible(False)
            continue

        worksheet = worksheets[axis_index]
        sheet_summary = well_summary[well_summary["worksheet"].eq(worksheet)]
        color = COLORS.get(worksheet, "#64748b")
        x = np.arange(sheet_summary.shape[0])
        if worksheet.endswith("-Lisa"):
            labels = [
                f"{iptg:g} uM\n{well}" if np.isfinite(float(iptg)) else str(well)
                for well, iptg in zip(sheet_summary["well"], sheet_summary["iptg_uM"])
            ]
            panel_title = f"{worksheet}: GFP BB BG by IPTG"
            x_label = "IPTG concentration / well"
        else:
            labels = sheet_summary["well"].astype(str).to_numpy()
            panel_title = f"{worksheet}: GFP BB BG by multiplex well"
            x_label = "Multiplex well/group"
        means = sheet_summary["mean"].to_numpy(dtype=float)
        sem = sheet_summary["sem"].fillna(0).to_numpy(dtype=float)

        ax.bar(
            x,
            means,
            yerr=sem,
            width=0.72,
            color=color,
            alpha=0.72,
            edgecolor="#1f2937",
            linewidth=0.8,
            capsize=3,
        )
        ax.set_title(panel_title)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylim(0, y_limit)
        ax.grid(True, axis="y", color="#d0d7de", linewidth=0.8, alpha=0.75)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel(x_label)
        if axis_index % 2 == 0:
            ax.set_ylabel("DonBB_BG mean +/- SEM")

    fig.suptitle(
        "GFP BB BG (DonBB_BG): LISA wells are IPTG titrations; Multi sheets are multiplexed",
        fontsize=16,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    svg_path = output_path.with_suffix(".svg")
    fig.savefig(svg_path)
    plt.close(fig)
    print(f"Wrote {output_path}")
    print(f"Wrote {svg_path}")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    workbook_path = prepare_workbook_path(args.workbook, output_dir, args.no_snapshot)
    print(f"Reading {workbook_path}")
    sheets = pd.read_excel(workbook_path, sheet_name=None, engine="openpyxl")

    raw_datasets: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    normalized_datasets: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    raw_results: dict[str, FitResult] = {}
    normalized_results: dict[str, FitResult] = {}
    well_summaries: list[pd.DataFrame] = []

    for worksheet, raw_df in sheets.items():
        validate_columns(worksheet, raw_df)
        well_summaries.append(summarize_gfp_by_well(worksheet, raw_df))
        df = finite_filtered_frame(raw_df.copy(), min_x=args.min_x)

        x_raw = df["AccBB_BG"].to_numpy(dtype=float)
        y_raw = df["CorrFRET"].to_numpy(dtype=float)
        donor_gfp = df["DonBB_BG"].to_numpy(dtype=float)
        y_norm = y_raw / donor_gfp
        norm_keep = (
            np.isfinite(x_raw)
            & np.isfinite(y_norm)
            & np.isfinite(donor_gfp)
            & (donor_gfp > 0)
            & (x_raw >= args.min_x)
        )

        raw_datasets[worksheet] = (x_raw, y_raw)
        normalized_datasets[worksheet] = (x_raw[norm_keep], y_norm[norm_keep])
        raw_results[worksheet] = fit_curve(
            "raw",
            worksheet,
            x_raw,
            y_raw,
            maxfev=args.maxfev,
        )
        normalized_results[worksheet] = fit_curve(
            "corrfret_over_donbb_bg",
            worksheet,
            x_raw[norm_keep],
            y_norm[norm_keep],
            maxfev=args.maxfev,
        )

    all_results = [*raw_results.values(), *normalized_results.values()]
    fit_table = pd.DataFrame([result_to_dict(result) for result in all_results])
    csv_path = output_dir / "fret_saturation_fit_parameters.csv"
    xlsx_path = output_dir / "fret_saturation_fit_parameters.xlsx"
    fit_table.to_csv(csv_path, index=False)
    fit_table.to_excel(xlsx_path, index=False)
    print(f"Wrote {csv_path}")
    print(f"Wrote {xlsx_path}")

    well_summary = pd.concat(well_summaries, ignore_index=True)
    well_summary_csv = output_dir / "gfp_bb_bg_by_well_summary.csv"
    well_summary_xlsx = output_dir / "gfp_bb_bg_by_well_summary.xlsx"
    well_summary.to_csv(well_summary_csv, index=False)
    well_summary.to_excel(well_summary_xlsx, index=False)
    print(f"Wrote {well_summary_csv}")
    print(f"Wrote {well_summary_xlsx}")

    plot_analysis(
        analysis="raw",
        datasets=raw_datasets,
        results=raw_results,
        output_path=output_dir / "fret_saturation_raw.png",
        x_label="Acceptor BB BG (AccBB_BG)",
        y_label="Corrected FRET (CorrFRET)",
        title="Raw FRET saturation fit, analysis_included = 1 and in_bleach_roi = 1",
        point_alpha=args.point_alpha,
        point_size=args.point_size,
    )

    plot_gfp_by_well(
        well_summary,
        output_path=output_dir / "gfp_bb_bg_by_well_barplots.png",
    )

    plot_analysis(
        analysis="corrfret_over_donbb_bg",
        datasets=normalized_datasets,
        results=normalized_results,
        output_path=output_dir / "fret_saturation_corrfret_over_donbb_bg.png",
        x_label="Acceptor BB BG (AccBB_BG)",
        y_label="CorrFRET / DonBB_BG",
        title="FRET response normalized by donor GFP, analysis_included = 1 and in_bleach_roi = 1",
        point_alpha=args.point_alpha,
        point_size=args.point_size,
    )


if __name__ == "__main__":
    main()
