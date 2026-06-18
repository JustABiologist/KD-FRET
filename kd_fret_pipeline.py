#!/usr/bin/env python3
"""
KD-FRET end-to-end processing pipeline.

This script consolidates the historical ImageJ macros (Plugins 1–3) and the
Python post-processing scripts (Plugins 4, 5, 8, 9) into a single CLI utility.
It orchestrates the following steps:

1. Registration + cropping of all measurement stacks via the ImageJ wrapper.
2. Extraction of a representative frame per FOV for Cellpose segmentation.
3. Optional execution of the Cellpose CLI to generate masks/ROIs.
4. Quantification of each registered stack using the Cellpose masks.
5. Computation of the donor/acceptor summary statistics used downstream.
6. Consolidation of the experiment into an appendable CSV plus an Excel workbook.

Assumptions that are in line with the original workflow:
* Measurement folders are named measurementX (case-insensitive) and contain
  xyNN subfolders or files tagged with xyNN (NN = 01..20).
* Each registered stack contains 54 frames (time points) saved as individual
  TIFFs (xyNNuse0001.tif ... xyNNuse0054.tif).
* Cellpose is invoked with --save_tif (so *_masks.tif exists) and the masks
  align with the registered stacks (cells do not move after registration).
* Donor/acceptor background values are user-provided for now (CLI arguments).
* FRET metrics reuse the same frame windows as Plugins 5/8/9.
* IPTG concentrations advance every N FOVs (default 20) across the entire dataset
  so downstream Kd fits can be automated. The last two FOVs of the FINAL
  measurement in the dataset are automatically labeled 'BLANK'.

Because the historical pipeline relied on several manual GUI steps, the script
prompts the user exactly once to draw the laser ROI inside ImageJ; the saved
ROI is reused for every subsequent stack.
"""

from __future__ import annotations

import argparse
import logging
import time
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops

try:
    import imagej
except ImportError as exc:  # pragma: no cover - handled at runtime
    raise ImportError(
        "pyimagej is required. Install the provided conda environment first."
    ) from exc


# --------------------------------------------------------------------------------------
# Dataclasses describing the measurement hierarchy
# --------------------------------------------------------------------------------------


@dataclass
class Measurement:
    """Metadata for a single measurement folder."""

    name: str
    source_dir: Path
    group: str
    index: int
    fovs: List[str] = field(default_factory=list)


@dataclass
class FovContext:
    """Concrete paths for each FOV across the pipeline stages."""

    measurement: Measurement
    label: str
    index: int
    iptg_concentration: Union[float, str]
    registered_stack_dir: Path
    cellpose_frame_path: Path
    cellpose_mask_path: Path


@dataclass
class PipelineConfig:
    input_root: Path
    output_root: Path
    registered_root: Path
    cellpose_input: Path
    results_csv: Path
    results_xlsx: Path
    backgrounds: Tuple[float, float]
    sequence_start: int
    cellpose_frame_index: int
    min_quality_ratio: float
    min_bleach_percent: float
    imagej_distribution: Optional[str]
    default_group: str
    reuse_laser_roi: Optional[Path]
    laser_roi_path: Path
    skip_registration: bool
    skip_cellpose: bool
    skip_measurement: bool
    cellpose_model: str
    cellpose_diameter: float
    cellpose_channels: Tuple[int, int]
    iptg_concentrations: List[float]
    fovs_per_iptg: int
    run_id: str
    csv_decimal: str


# --------------------------------------------------------------------------------------
# CLI parsing helpers
# --------------------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified KD-FRET processing pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Folder containing measurement sub-folders (measurement1, ...).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("fret_output"),
        help="Root folder for all intermediate + final artefacts.",
    )
    parser.add_argument(
        "--background-donor",
        type=float,
        required=True,
        help="Background value for the donor channel (BG don).",
    )
    parser.add_argument(
        "--background-acceptor",
        type=float,
        required=True,
        help="Background value for the acceptor channel (BG acc).",
    )
    parser.add_argument(
        "--sequence-start",
        type=int,
        default=5,
        help="Start index used by File.openSequence (same as Plugin 1).",
    )
    parser.add_argument(
        "--cellpose-frame-index",
        type=int,
        default=5,
        help="Frame number extracted per FOV for Cellpose (1-based).",
    )
    parser.add_argument(
        "--min-quality-ratio",
        type=float,
        default=0.8,
        help="Minimum Mean/StdDev ratio to accept a ROI.",
    )
    parser.add_argument(
        "--min-bleach-percent",
        type=float,
        default=69.6,
        help="Bleach percentage threshold used in historical DaySummaryIN.",
    )
    parser.add_argument(
        "--imagej-distribution",
        type=str,
        default=None,
        help="Optional Fiji/ImageJ distribution string passed to imagej.init().",
    )
    parser.add_argument(
        "--default-group",
        type=str,
        default="IN",
        help="Fallback label if a measurement name does not encode IN/OUT.",
    )
    parser.add_argument(
        "--laser-roi",
        type=Path,
        default=None,
        help="Reuse an existing ROI instead of prompting inside ImageJ.",
    )
    parser.add_argument(
        "--skip-registration",
        action="store_true",
        help="Assume registered stacks already exist.",
    )
    parser.add_argument(
        "--skip-cellpose",
        action="store_true",
        help="Do not invoke Cellpose (expects *_masks.tif to exist).",
    )
    parser.add_argument(
        "--skip-measurement",
        action="store_true",
        help="Skip quantitative analysis (useful when only preparing data).",
    )
    parser.add_argument(
        "--cellpose-model",
        type=str,
        default="bact_fluor_cp3",
        help="Cellpose pretrained model identifier.",
    )
    parser.add_argument(
        "--cellpose-diameter",
        type=float,
        default=25.0,
        help="Expected object diameter passed to Cellpose.",
    )
    parser.add_argument(
        "--cellpose-channels",
        type=int,
        nargs=2,
        default=(0, 0),
        metavar=("CHAN", "CHAN2"),
        help="Cellpose channel specification, e.g. 0 0 for grayscale.",
    )
    parser.add_argument(
        "--existing-results",
        type=Path,
        default=None,
        help="Path to an existing CSV that should be appended.",
    )
    parser.add_argument(
        "--results-xlsx",
        type=Path,
        default=None,
        help="Optional Excel file path; defaults to <output_root>/fret_results.xlsx.",
    )
    parser.add_argument(
        "--iptg-concentrations",
        type=float,
        nargs="+",
        required=True,
        help="List of IPTG concentrations consumed sequentially across all FOVs.",
    )
    parser.add_argument(
        "--fovs-per-iptg",
        type=int,
        default=20,
        help="Number of FOVs per concentration block (global order).",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        default=datetime.now().strftime("%Y%m%d_%H%M%S"),
        help="Run identifier stored in the final CSV.",
    )
    parser.add_argument(
        "--csv-decimal",
        type=str,
        choices=[".", ","],
        default=".",
        help="Decimal separator to use when reading/writing CSV outputs.",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Verbosity for console logging.",
    )
    return parser.parse_args()


# --------------------------------------------------------------------------------------
# Utility helpers
# --------------------------------------------------------------------------------------


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)-8s | %(message)s",
    )


def infer_group(name: str, default_group: str) -> str:
    lowered = name.lower()
    if "out" in lowered:
        return "OUT"
    if "in" in lowered:
        return "IN"
    return default_group.upper()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def iptg_for_fov(
    *,
    local_index: int,
    global_index: int,
    concentrations: Sequence[float],
    block_size: int,
    total_fovs_in_measurement: int,
    is_last_measurement: bool,
) -> Union[float, str]:
    """Determine IPTG value for a FOV.

    IPTG concentrations advance sequentially across the entire dataset (global_index),
    but the BLANK logic still depends on the final measurement (local_index).
    """
    # The last two FOVs of the LAST measurement are reserved as BLANK
    if is_last_measurement and local_index > total_fovs_in_measurement - 2:
        return "BLANK"

    block = (global_index - 1) // block_size
    if block < 0 or block >= len(concentrations):
        raise ValueError(
            f"Global FOV index {global_index} exceeds provided IPTG blocks. "
            f"Blocks available: {len(concentrations)} with size {block_size}."
        )
    return concentrations[block]


def path_for_macro(path: Path) -> str:
    """Convert a filesystem path into an ImageJ-friendly string."""
    return str(path).replace("\\", "/")


def detect_fovs(measurement_dir: Path) -> List[str]:
    """Collect xyNN labels present inside a measurement folder."""
    fov_dirs = [
        p.name
        for p in measurement_dir.iterdir()
        if p.is_dir() and re.match(r"xy\d{2}", p.name.lower())
    ]
    if fov_dirs:
        return sorted(fov_dirs)

    pattern = re.compile(r"(xy\d{2})", re.IGNORECASE)
    labels: set[str] = set()
    for file in measurement_dir.glob("**/*.tif"):
        match = pattern.search(file.name)
        if match:
            labels.add(match.group(1))
    return sorted(labels)


def collect_measurements(input_root: Path, default_group: str) -> List[Measurement]:
    measurements: List[Measurement] = []
    for idx, folder in enumerate(sorted(p for p in input_root.iterdir() if p.is_dir()), start=1):
        labels = detect_fovs(folder)
        if not labels:
            logging.warning("Skipping %s because no xyNN FOVs were detected.", folder)
            continue
        name = f"Me{idx:02d}"
        measurement = Measurement(
            name=name,
            source_dir=folder,
            group=infer_group(folder.name, default_group),
            index=idx,
            fovs=labels,
        )
        measurements.append(measurement)
    if not measurements:
        raise ValueError(f"No measurement folders found inside {input_root}.")
    return measurements


# --------------------------------------------------------------------------------------
# ImageJ integration
# --------------------------------------------------------------------------------------


SIFT_PARAMS = (
    "initial_gaussian_blur=1.60 steps_per_scale_octave=3 "
    "minimum_image_size=64 maximum_image_size=512 "
    "feature_descriptor_size=4 feature_descriptor_orientation_bins=8 "
    "closest/next_closest_ratio=0.85 maximal_alignment_error=10 "
    "inlier_ratio=0.10 expected_transformation=Rigid interpolate "
    "show_info"
)


def init_imagej(distribution: Optional[str]) -> imagej.ImageJ:
    logging.info("Starting ImageJ (distribution=%s)...", distribution or "sc.fiji:fiji")
    ij = imagej.init(distribution or "sc.fiji:fiji", headless=False)
    return ij


def align_stack_only(
    ij: imagej.ImageJ,
    measurement: Measurement,
    config: PipelineConfig,
    fov: str,
    dest_dir: Path,
) -> Path:
    """
    Align a single stack without cropping and persist the registered frames at dest_dir.
    """
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    dest_dir = ensure_dir(dest_dir)
    logging.info(
        "Preparing ROI sample by aligning %s %s into %s",
        measurement.name,
        fov,
        dest_dir,
    )
    macro = f"""
run("Close All");
src = "{path_for_macro(measurement.source_dir)}";
File.openSequence(src, "filter={fov} start={config.sequence_start}");
run("Enhance Contrast", "saturated=0.35");
print("Loaded stack has " + nSlices + " slices, starting SIFT alignment...");
run("Linear Stack Alignment with SIFT", "{SIFT_PARAMS}");
wait(4000);
alignedTitle = "Aligned 54 of 54";
if (isOpen(alignedTitle)) {{
    selectWindow(alignedTitle);
}} else if (isOpen("Aligned")) {{
    selectWindow("Aligned");
}} else {{
    print("ERROR: Aligned stack window 'Aligned 54 of 54' not found after SIFT");
    exit();
}}
print("Aligned stack active: " + getTitle());
dest = "{path_for_macro(dest_dir)}";
if (!File.exists(dest)) {{
    File.makeDirectory(dest);
}}
print("Saving aligned sequence to " + dest + "/");
run("Image Sequence... ", "select=[" + dest + "/] dir=[" + dest + "/] format=TIFF name={fov}use");
print("Saved aligned sequence to " + dest + "/");
list = getFileList(dest + "/");
print("Verified " + list.length + " files in destination.");
run("Close All");
"""
    ij.py.run_macro(macro)
    tiffs = list(dest_dir.glob("*.tif"))
    logging.debug("Aligned stack files: %s", [f.name for f in tiffs])
    if not tiffs:
        raise RuntimeError(
            f"Failed to materialize aligned stack for {measurement.name} {fov} at {dest_dir}"
        )
    logging.info("Aligned stack materialized: %s (%d TIFFs)", dest_dir, len(tiffs))
    return dest_dir


def prompt_for_laser_roi_from_stack(
    ij: imagej.ImageJ,
    stack_dir: Path,
    config: PipelineConfig,
) -> None:
    """
    Open an already aligned stack (mCherry channel) and let the user draw the bleaching ROI.
    """
    logging.info("Opening aligned stack at %s for ROI definition.", stack_dir)
    macro = f"""
run("Close All");
src = "{path_for_macro(stack_dir)}";
File.openSequence(src, "start=1");
setSlice(nSlices);
run("Enhance Contrast", "saturated=0.40");
run("Brightness/Contrast...");
setTool("oval");
waitForUser("Laser ROI", "Draw the bleaching ROI on this aligned stack (Frame " + nSlices + ").\\n\\nUse the Toolbar to switch to Rectangle if needed.\\nAdjust Brightness/Contrast if needed.\\nThen press OK.");
roiManager("Reset");
roiManager("Add");
roiManager("Select", 0);
roiManager("Save", "{path_for_macro(config.laser_roi_path)}");
run("Close All");
"""
    ij.py.run_macro(macro)
    logging.info("Saved laser ROI to %s", config.laser_roi_path)


def register_and_crop(
    ij: imagej.ImageJ,
    measurement: Measurement,
    config: PipelineConfig,
    fov: str,
) -> Path:
    dest_parent = ensure_dir(config.registered_root / measurement.name)
    dest_dir = ensure_dir(dest_parent / fov)
    macro = f"""
run("Close All");
src = "{path_for_macro(measurement.source_dir)}";
File.openSequence(src, "filter={fov} start={config.sequence_start}");
run("Enhance Contrast", "saturated=0.35");
print("Loaded stack has " + nSlices + " slices, starting SIFT alignment...");
run("Linear Stack Alignment with SIFT", "{SIFT_PARAMS}");
wait(4000);
alignedTitle = "Aligned 54 of 54";
if (isOpen(alignedTitle)) {{
    selectWindow(alignedTitle);
}} else if (isOpen("Aligned")) {{
    selectWindow("Aligned");
}} else {{
    print("ERROR: Aligned stack window 'Aligned 54 of 54' not found after SIFT");
    exit();
}}
print("Aligned stack active: " + getTitle());
roiManager("Reset");
roiManager("Open", "{path_for_macro(config.laser_roi_path)}");
roiManager("Select", 0);
run("Crop");
run("Enhance Contrast", "saturated=0.35");
resetMinAndMax();
    dest = "{path_for_macro(dest_dir)}";
    if (!File.exists(dest)) {{
        File.makeDirectory(dest);
    }}
    print("Saving cropped sequence to " + dest + "/");
    run("Image Sequence... ", "select=[" + dest + "/] dir=[" + dest + "/] format=TIFF name={fov}use");
    print("Saved cropped sequence to " + dest + "/");
    list = getFileList(dest + "/");
    print("Verified " + list.length + " files in destination.");
run("Close All");
"""
    ij.py.run_macro(macro)
    return dest_dir


# --------------------------------------------------------------------------------------
# Cellpose helpers
# --------------------------------------------------------------------------------------


def extract_cellpose_frames(
    measurement: Measurement,
    fov: str,
    stack_dir: Path,
    cellpose_dir: Path,
    frame_index: int,
) -> Tuple[Path, Path]:
    files = sorted(stack_dir.glob("*.tif"))
    if not files:
        raise FileNotFoundError(f"No TIFFs found in {stack_dir}")
    if frame_index < 1 or frame_index > len(files):
        raise ValueError(
            f"Frame index {frame_index} outside [1, {len(files)}] for {stack_dir}"
        )
    selected = files[frame_index - 1]
    target_name = f"{measurement.name}{fov}use{frame_index:04d}.tif"
    dest = cellpose_dir / target_name
    shutil.copyfile(selected, dest)
    mask_path = dest.with_name(dest.stem + "_cp_masks.tif")
    return dest, mask_path


def run_cellpose_cli(
    cellpose_input: Path,
    model: str,
    diameter: float,
    channels: Tuple[int, int],
) -> None:
    cmd = [
        sys.executable,
        "-m",
        "cellpose",
        "--dir",
        str(cellpose_input),
        "--pretrained_model",
        model,
        "--diameter",
        str(diameter),
        "--save_tif",
        "--save_png",
        "--save_outlines",
        "--save_rois",
        "--verbose",
        "--no_npy",
        "--chan",
        str(channels[0]),
        "--chan2",
        str(channels[1]),
    ]
    logging.info("Running Cellpose: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


# --------------------------------------------------------------------------------------
# Quantification helpers
# --------------------------------------------------------------------------------------


def load_stack(stack_dir: Path) -> np.ndarray:
    frames = sorted(stack_dir.glob("*.tif"))
    if not frames:
        raise FileNotFoundError(f"No TIFF files inside {stack_dir}")
    stack = np.stack([tifffile.imread(frame) for frame in frames], axis=0)
    return stack


def load_mask(mask_path: Path) -> np.ndarray:
    if not mask_path.exists():
        raise FileNotFoundError(mask_path)
    mask = tifffile.imread(mask_path)
    if mask.ndim == 3:
        mask = mask[0]
    return mask


def imagej_std(values: np.ndarray, axis: int) -> np.ndarray:
    """Match ImageJ Results-table StdDev: sample stddev, zero for one-pixel ROIs."""
    if values.shape[axis] <= 1:
        out_shape = list(values.shape)
        out_shape.pop(axis)
        return np.zeros(out_shape, dtype=float)
    return np.std(values, axis=axis, ddof=1)


def compute_roi_timeseries(
    stack: np.ndarray,
    mask: np.ndarray,
    measurement: Measurement,
    fov: str,
    stack_dir: Path,
    mask_path: Path,
    fov_index: int,
    iptg_concentration: Union[float, str],
    config: PipelineConfig,
) -> List[Dict[str, object]]:
    if stack.shape[1:] != mask.shape:
        raise ValueError(
            f"Mask mismatch for {measurement.name}-{fov}: stack {stack.shape[1:]}, mask {mask.shape}"
        )
    props = {prop.label: prop for prop in regionprops(mask)}
    labels = sorted(label for label in np.unique(mask) if label > 0)
    records: List[Dict[str, object]] = []
    bg_don, bg_acc = config.backgrounds
    frame_numbers = np.arange(1, stack.shape[0] + 1)
    for label in labels:
        cell_mask = mask == label
        area_px = int(cell_mask.sum())
        prop = props.get(label)
        if prop is None:
            continue
        angle = float(np.degrees(prop.orientation))
        major = float(prop.major_axis_length)
        minor = float(prop.minor_axis_length)
        series = stack[:, cell_mask]
        mean_series = series.mean(axis=1)
        std_series = imagej_std(series, axis=1)

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_series = np.where(std_series != 0, mean_series / std_series, np.nan)
        avg_ratio = np.nanmean(ratio_series)
        if np.isnan(avg_ratio) or avg_ratio < config.min_quality_ratio:
            logging.debug(
                "Discarding %s-%s label %s due to quality ratio %.3f",
                measurement.name,
                fov,
                label,
                avg_ratio,
            )
            continue

        def avg_range(start: int, end: int) -> float:
            sl = slice(start - 1, end)
            return float(np.nanmean(mean_series[sl]))

        # Corrected Frame Indices (1-based) based on Legacy Plugin logic:
        # DonBB (Baseline): Frames 23-25 (was 26-28)
        # DonAB (Signal):   Frames 26-28 (was 29-31)
        # AccBB (Baseline): Frames 51-52 (was 2-3)
        # AccAB (Bleached): Frames 53-54 (was 54-55)
        don_bb = avg_range(23, 25)
        don_ab = avg_range(26, 28)
        acc_bb = avg_range(51, 52)
        acc_ab = avg_range(53, 54)

        don_bb_bg = don_bb - bg_don
        don_ab_bg = don_ab - bg_don
        acc_bb_bg = acc_bb - bg_acc
        acc_ab_bg = acc_ab - bg_acc

        fret = np.nan
        corr_fret = np.nan
        bleached = np.nan
        fac = np.nan
        ratio_metric = np.nan

        if don_ab_bg != 0:
            fret = 100.0 * (don_ab_bg - don_bb_bg) / don_ab_bg

        if acc_bb_bg != 0:
            bleached = 100.0 * (acc_bb_bg - acc_ab_bg) / acc_bb_bg

        if not np.isnan(bleached):
            fac = (100.0 - bleached) / 100.0

        denominator = don_ab_bg - (fac * don_bb_bg if not np.isnan(fac) else 0.0)
        if denominator != 0:
            corr_fret = 100.0 * (don_ab_bg - don_bb_bg) / denominator

        if don_ab_bg != 0:
            ratio_metric = acc_bb_bg / don_ab_bg

        record = {
            "run_id": config.run_id,
            "measurement": measurement.name,
            "measurement_group": measurement.group,
            "fov": fov,
            "fov_index": fov_index,
            "iptg_concentration": iptg_concentration,
            "cell_label": int(label),
            "frame_count": int(stack.shape[0]),
            "area_px": area_px,
            "major_axis_px": major,
            "minor_axis_px": minor,
            "angle_deg": angle,
            "quality_ratio": float(avg_ratio),
            "DonBB": don_bb,
            "DonAB": don_ab,
            "AccBB": acc_bb,
            "AccAB": acc_ab,
            "BG_don": bg_don,
            "BG_acc": bg_acc,
            "DonBB_BG": don_bb_bg,
            "DonAB_BG": don_ab_bg,
            "AccBB_BG": acc_bb_bg,
            "ACCAB_BG": acc_ab_bg,
            "FRET": fret,
            "BleachedPercent": bleached,
            "Fac": fac,
            "CorrFRET": corr_fret,
            "Ratio": ratio_metric,
            "BleachedPass": (
                False
                if np.isnan(bleached)
                else bool(bleached >= config.min_bleach_percent)
            ),
            "source_stack_dir": str(stack_dir),
            "cellpose_mask_path": str(mask_path),
        }
        records.append(record)
    return records


def process_measurements(
    contexts: List[FovContext],
    config: PipelineConfig,
) -> pd.DataFrame:
    all_records: List[Dict[str, object]] = []
    for ctx in contexts:
        if not ctx.registered_stack_dir.exists():
            logging.warning("Missing registered stack for %s %s", ctx.measurement.name, ctx.label)
            continue
        if not ctx.cellpose_mask_path.exists():
            logging.warning("Missing Cellpose mask for %s %s", ctx.measurement.name, ctx.label)
            continue
        stack = load_stack(ctx.registered_stack_dir)
        mask = load_mask(ctx.cellpose_mask_path)
        records = compute_roi_timeseries(
            stack=stack,
            mask=mask,
            measurement=ctx.measurement,
            fov=ctx.label,
            stack_dir=ctx.registered_stack_dir,
            mask_path=ctx.cellpose_mask_path,
            fov_index=ctx.index,
            iptg_concentration=ctx.iptg_concentration,
            config=config,
        )
        all_records.extend(records)
    if not all_records:
        raise RuntimeError("No valid ROIs were quantified. Check Cellpose masks/quality filters.")
    df = pd.DataFrame(all_records)
    return df


def append_and_save(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    ensure_dir(config.results_csv.parent)
    if config.results_csv.exists():
        logging.info("Appending to existing results at %s", config.results_csv)
        previous = pd.read_csv(config.results_csv, decimal=config.csv_decimal)
        df = pd.concat([previous, df], ignore_index=True)
    df.to_csv(config.results_csv, index=False, decimal=config.csv_decimal)
    logging.info("Saved %s rows to %s", len(df), config.results_csv)
    return df


def save_excel(df: pd.DataFrame, config: PipelineConfig) -> Path:
    ensure_dir(config.results_xlsx.parent)
    with pd.ExcelWriter(config.results_xlsx, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Results")
    logging.info("Saved %s rows to %s", len(df), config.results_xlsx)
    return config.results_xlsx


# --------------------------------------------------------------------------------------
# Main orchestration
# --------------------------------------------------------------------------------------


def build_fov_contexts(
    measurements: List[Measurement],
    config: PipelineConfig,
) -> List[FovContext]:
    contexts: List[FovContext] = []
    global_fov_index = 1
    num_measurements = len(measurements)
    for m_idx, measurement in enumerate(measurements):
        is_last_measurement = (m_idx == num_measurements - 1)
        for idx, fov in enumerate(measurement.fovs, start=1):
            current_global_index = global_fov_index
            global_fov_index += 1
            registered_dir = config.registered_root / measurement.name / fov
            try:
                cellpose_frame, mask_path = extract_cellpose_frames(
                    measurement=measurement,
                    fov=fov,
                    stack_dir=registered_dir,
                    cellpose_dir=config.cellpose_input,
                    frame_index=config.cellpose_frame_index,
                )
            except (FileNotFoundError, ValueError) as exc:
                logging.warning(
                    "Skipping %s %s because the registered stack is missing: %s",
                    measurement.name,
                    fov,
                    exc,
                )
                continue
            try:
                iptg_value = iptg_for_fov(
                    local_index=idx,
                    global_index=current_global_index,
                    concentrations=config.iptg_concentrations,
                    block_size=config.fovs_per_iptg,
                    total_fovs_in_measurement=len(measurement.fovs),
                    is_last_measurement=is_last_measurement,
                )
            except ValueError as exc:
                raise ValueError(
                    f"Measurement {measurement.name}: {exc}"
                ) from exc
            ctx = FovContext(
                measurement=measurement,
                label=fov,
                index=idx,
                iptg_concentration=iptg_value,
                registered_stack_dir=registered_dir,
                cellpose_frame_path=cellpose_frame,
                cellpose_mask_path=mask_path,
            )
            contexts.append(ctx)
    return contexts


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    config = PipelineConfig(
        input_root=args.input_root.resolve(),
        output_root=args.output_root.resolve(),
        registered_root=(args.output_root / "01_registered").resolve(),
        cellpose_input=(args.output_root / "02_cellpose_input").resolve(),
        results_csv=(
            args.existing_results.resolve()
            if args.existing_results
            else (args.output_root / "fret_results.csv").resolve()
        ),
        results_xlsx=(
            args.results_xlsx.resolve()
            if args.results_xlsx
            else (args.output_root / "fret_results.xlsx").resolve()
        ),
        backgrounds=(args.background_donor, args.background_acceptor),
        sequence_start=args.sequence_start,
        cellpose_frame_index=args.cellpose_frame_index,
        min_quality_ratio=args.min_quality_ratio,
        min_bleach_percent=args.min_bleach_percent,
        imagej_distribution=args.imagej_distribution,
        default_group=args.default_group,
        reuse_laser_roi=args.laser_roi,
        laser_roi_path=(args.laser_roi or (args.output_root / "laser_roi.roi")).resolve(),
        skip_registration=args.skip_registration,
        skip_cellpose=args.skip_cellpose,
        skip_measurement=args.skip_measurement,
        cellpose_model=args.cellpose_model,
        cellpose_diameter=args.cellpose_diameter,
        cellpose_channels=tuple(args.cellpose_channels),
        iptg_concentrations=args.iptg_concentrations,
        fovs_per_iptg=args.fovs_per_iptg,
        run_id=args.run_id,
        csv_decimal=args.csv_decimal,
    )

    ensure_dir(config.output_root)
    ensure_dir(config.registered_root)
    ensure_dir(config.cellpose_input)

    measurements = collect_measurements(config.input_root, config.default_group)
    logging.info("Discovered %d measurements.", len(measurements))

    ij = None
    if not config.skip_registration or config.reuse_laser_roi is None:
        ij = init_imagej(config.imagej_distribution)

    if config.reuse_laser_roi and config.reuse_laser_roi.exists():
        config.laser_roi_path = config.reuse_laser_roi.resolve()
        logging.info("Using existing ROI at %s", config.laser_roi_path)
    elif not config.skip_registration:
        roi_sample_dir = config.registered_root / "roi_temp"
        sample_found = False
        for sample_measurement in measurements:
            for sample_fov in sample_measurement.fovs:
                try:
                    align_stack_only(
                        ij=ij,
                        measurement=sample_measurement,
                        config=config,
                        fov=sample_fov,
                        dest_dir=roi_sample_dir,
                    )
                    prompt_for_laser_roi_from_stack(ij, roi_sample_dir, config)
                    sample_found = True
                    break
                except RuntimeError as exc:
                    logging.warning(
                        "Skipping %s %s for ROI setup: %s",
                        sample_measurement.name,
                        sample_fov,
                        exc,
                    )
                    continue
                finally:
                    if roi_sample_dir.exists():
                        shutil.rmtree(roi_sample_dir, ignore_errors=True)
            if sample_found:
                break
        if not sample_found:
            raise RuntimeError(
                "Unable to align any FOV for ROI definition; check image contrast or SIFT settings."
            )
    else:
        raise ValueError("Cannot skip registration without providing --laser-roi.")

    if not config.skip_registration:
        for measurement in measurements:
            for fov in measurement.fovs:
                logging.info("Registering %s %s", measurement.name, fov)
                register_and_crop(ij, measurement, config, fov)
            logging.info("Pausing 10 seconds after %s to mimic Plugin 1 behavior.", measurement.name)
            time.sleep(10)
    else:
        logging.info("Skipping registration as requested.")

    contexts = build_fov_contexts(measurements, config)
    logging.info("Prepared %d FOV contexts for Cellpose.", len(contexts))

    if not config.skip_cellpose:
        run_cellpose_cli(
            cellpose_input=config.cellpose_input,
            model=config.cellpose_model,
            diameter=config.cellpose_diameter,
            channels=config.cellpose_channels,
        )
    else:
        logging.info("Skipping Cellpose execution.")

    if config.skip_measurement:
        logging.info("Measurement step skipped; pipeline finished early.")
        return

    df = process_measurements(contexts, config)
    combined_df = append_and_save(df, config)
    save_excel(combined_df, config)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover
        logging.error("Pipeline failed: %s", exc)
        sys.exit(1)

