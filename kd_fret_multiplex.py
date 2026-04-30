#!/usr/bin/env python3
"""
KD-FRET multiplex processing pipeline.

This script keeps the legacy registered-TIFF multiplex workflow, and adds a raw
ND2 workflow for experiments where each measurement is a subdirectory containing
files named like:

    YYYYMMDD_HHMMSS_001_<channelname>_Seq0000.nd2

The raw workflow maps Seq0000..Seq0004 to the acquisition steps, groups files by
measurement folder and ND2 position label,
creates Cellpose inputs, quantifies cells in the bleached mask and the 100 px
surrounding unbleached ring, and writes one workbook with one sheet per
measurement folder.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops
from skimage import draw, filters, measure, morphology
from skimage.restoration import rolling_ball

try:
    import imagej
except ImportError as exc:  # pragma: no cover - handled at runtime
    imagej = None


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
    registered_bleached_root: Path
    registered_unbleached_root: Path
    cellpose_input: Path
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
    iptg_concentrations: Optional[List[float]]
    fovs_per_iptg: int
    run_id: str
    csv_decimal: str
    multiplex: bool


@dataclass(frozen=True)
class RawNd2File:
    path: Path
    timestamp: str
    token: str
    channel: str
    role: str
    seq: int


@dataclass
class RawFovJob:
    measurement_name: str
    measurement_dir: Path
    position_index: int
    fov: str
    well_index: int
    well: str
    iptg_label: str
    condition: str
    role_files: Dict[str, RawNd2File]
    cellpose_frame_path: Path
    cellpose_mask_path: Path


# --------------------------------------------------------------------------------------
# CLI parsing helpers
# --------------------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified KD-FRET Multiplex processing pipeline.",
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
        default=Path("fret_output_multiplex"),
        help="Root folder for all intermediate + final artefacts.",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "raw-nd2", "legacy-tiff"],
        default="auto",
        help="Use raw ND2 multiplex ingestion, legacy registered TIFF ingestion, or auto-detect.",
    )
    parser.add_argument(
        "--output-xlsx",
        type=Path,
        default=None,
        help="Workbook path for raw ND2 multiplex output. Defaults to <output_root>/multiplex_results.xlsx.",
    )
    parser.add_argument(
        "--background-donor",
        type=float,
        default=None,
        help="Background value for the donor channel (BG don).",
    )
    parser.add_argument(
        "--background-acceptor",
        type=float,
        default=None,
        help="Background value for the acceptor channel (BG acc).",
    )
    parser.add_argument(
        "--background-mode",
        choices=["manual", "auto", "rolling-ball"],
        default="auto",
        help=(
            "Background correction for raw ND2 mode. manual uses --background-donor/acceptor; "
            "auto estimates scalar background from cell-free pixels; rolling-ball subtracts "
            "a per-frame rolling-ball background first."
        ),
    )
    parser.add_argument(
        "--rolling-ball-radius",
        type=int,
        default=50,
        help="Rolling-ball radius in pixels for --background-mode rolling-ball.",
    )
    parser.add_argument(
        "--background-exclusion-radius",
        type=int,
        default=8,
        help="Dilate cell masks by this many pixels before estimating cell-free background.",
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
        "--donor-channel-name",
        type=str,
        default=None,
        help="Substring used to identify donor-channel ND2 files. If omitted, aliases are inferred.",
    )
    parser.add_argument(
        "--acceptor-channel-name",
        type=str,
        default=None,
        help="Substring used to identify acceptor-channel ND2 files. If omitted, aliases are inferred.",
    )
    parser.add_argument(
        "--laser-channel-name",
        type=str,
        default=None,
        help="Substring used to identify laser/bleach ND2 files. If omitted, aliases are inferred.",
    )
    parser.add_argument(
        "--role-order",
        nargs=5,
        default=["donor_before", "acceptor_before", "laser", "acceptor_after", "donor_after"],
        choices=["donor_before", "donor_after", "acceptor_before", "acceptor_after", "laser"],
        help=(
            "Fallback acquisition role order for the five ND2 files when Seq numbers are not 0-4."
        ),
    )
    parser.add_argument(
        "--fov-labels",
        nargs="+",
        default=None,
        help="FOV/position labels in ND2 position order. These replace filename-based FOV inference.",
    )
    parser.add_argument(
        "--measurement-labels-json",
        type=str,
        default=None,
        help=(
            "JSON string or JSON file path with per-measurement well/condition labels. "
            "Example: '[[\"0uM\", \"3uM\", \"10uM\"], [\"A\", \"B\"]]'. "
            "Sublists follow sorted measurement subdirectory order."
        ),
    )
    parser.add_argument(
        "--fovs-per-well-by-measurement",
        type=int,
        nargs="+",
        default=None,
        help=(
            "FOVs per well for each sorted measurement subdirectory. "
            "Example: 5 10 10 5."
        ),
    )
    parser.add_argument(
        "--donor-before-count",
        type=int,
        default=2,
        help="Raw ND2 donor frames before bleaching, taken from the first donor frames per FOV.",
    )
    parser.add_argument(
        "--donor-after-count",
        type=int,
        default=2,
        help="Raw ND2 donor frames after bleaching, taken from the last donor frames per FOV.",
    )
    parser.add_argument(
        "--acceptor-before-count",
        type=int,
        default=20,
        help="Raw ND2 acceptor frames before bleaching, taken from the first acceptor frames per FOV.",
    )
    parser.add_argument(
        "--acceptor-after-count",
        type=int,
        default=18,
        help="Raw ND2 acceptor frames after bleaching, taken from the last acceptor frames per FOV.",
    )
    parser.add_argument(
        "--laser-count",
        type=int,
        default=4,
        help="Raw ND2 laser/bleach frames per FOV, used for automatic bleach-mask detection.",
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
        default=70.0,
        help="Bleach percentage threshold used in historical DaySummaryIN.",
    )
    parser.add_argument(
        "--roi-buffer-px",
        type=int,
        default=100,
        help="Pixels around the bleached mask to classify as the unbleached surrounding area.",
    )
    parser.add_argument(
        "--roi-mode",
        choices=["prompt-ring", "prompt-two", "auto-laser"],
        default="prompt-ring",
        help=(
            "Raw ND2 ROI mode: prompt bleached ROI and generate a 100 px ring, prompt bleached "
            "and unbleached ROIs, or infer the bleached ROI from the laser stack."
        ),
    )
    parser.add_argument(
        "--min-area-overlap",
        type=float,
        default=0.5,
        help="Minimum fraction of a cell mask that must overlap bleach/ring masks for classification.",
    )
    parser.add_argument(
        "--laser-threshold-percentile",
        type=float,
        default=99.5,
        help="Percentile threshold used when deriving a bleach mask from the laser channel.",
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
        help="Reuse an existing ImageJ ROI instead of prompting/detecting the laser area.",
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
    # In Multiplex mode, this is optional/ignored
    parser.add_argument(
        "--iptg-concentrations",
        type=float,
        nargs="+",
        default=None,
        help="List of IPTG concentrations (Ignored if --multiplex is active).",
    )
    parser.add_argument(
        "--fov-map",
        type=Path,
        default=None,
        help="Optional CSV with fov, well, iptg_label/iptg, and/or condition/name columns.",
    )
    parser.add_argument(
        "--fovs-per-well",
        type=int,
        default=0,
        help="Single FOVs-per-well value reused for every measurement in raw ND2 mode.",
    )
    parser.add_argument(
        "--fovs-per-iptg",
        type=int,
        default=20,
        help="Number of FOVs per concentration block (Ignored if --multiplex is active).",
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
    parser.add_argument(
        "--multiplex",
        action="store_true",
        help="Enable multiplex mode: Output per measurement, ignore IPTG blocks.",
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


def path_for_macro(path: Path) -> str:
    """Convert a filesystem path into an ImageJ-friendly string."""
    return str(path).replace("\\", "/")


FRAME_FILE_RE = re.compile(r"use(\d{4})\.tiff?$", re.IGNORECASE)
ND2_FILE_RE = re.compile(
    r"(?P<timestamp>\d{6,8}_\d{6})_(?P<token>\d{3})_(?P<channel>.+)_Seq(?P<seq>\d+)\.nd2$",
    re.IGNORECASE,
)


def stack_frame_files(stack_dir: Path) -> List[Path]:
    """Return image-sequence frames only, excluding mask/ROI sidecar TIFFs."""
    indexed: List[Tuple[int, Path]] = []
    for path in stack_dir.glob("*.tif*"):
        match = FRAME_FILE_RE.search(path.name)
        if match:
            indexed.append((int(match.group(1)), path))
    if indexed:
        return [path for _, path in sorted(indexed)]

    sidecar_terms = ("mask", "roi", "outline")
    return sorted(
        path
        for path in stack_dir.glob("*.tif*")
        if not any(term in path.stem.lower() for term in sidecar_terms)
    )


def has_raw_nd2(input_root: Path) -> bool:
    return any(input_root.rglob("*.nd2"))


def should_use_raw_nd2(args: argparse.Namespace) -> bool:
    if args.mode == "raw-nd2":
        return True
    if args.mode == "legacy-tiff":
        return False
    return has_raw_nd2(args.input_root)


def sanitize_sheet_name(name: str, used: set[str]) -> str:
    cleaned = re.sub(r"[\[\]\:\*\?\/\\]", "_", name).strip() or "Measurement"
    cleaned = cleaned[:31]
    candidate = cleaned
    suffix = 1
    while candidate in used:
        tail = f"_{suffix}"
        candidate = cleaned[: 31 - len(tail)] + tail
        suffix += 1
    used.add(candidate)
    return candidate


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
# Raw ND2 multiplex helpers
# --------------------------------------------------------------------------------------


CHANNEL_ALIASES: Mapping[str, Tuple[str, ...]] = {
    "donor_before": ("donorbefore", "donorpre", "donbefore", "donpre", "gfpbefore", "gfppre"),
    "donor_after": ("donorafter", "donorpost", "donafter", "donpost", "gfpafter", "gfppost"),
    "acceptor_before": ("acceptorbefore", "acceptorpre", "accbefore", "accpre", "mcherrybefore", "mcherrypre"),
    "acceptor_after": ("acceptorafter", "acceptorpost", "accafter", "accpost", "mcherryafter", "mcherrypost"),
    "donor": ("donor", "don", "gfp", "cfp", "fitc", "green"),
    "acceptor": ("acceptor", "acc", "mcherry", "rfp", "tritc", "cy3", "cy5", "red"),
    "laser": ("laser", "bleach", "acal", "photo", "405"),
}

SEQ_ROLE_MAP: Mapping[int, str] = {
    0: "donor_before",
    1: "acceptor_before",
    2: "laser",
    3: "acceptor_after",
    4: "donor_after",
}


def normalize_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def infer_channel_role(channel_name: str, args: argparse.Namespace) -> Optional[str]:
    normalized = normalize_token(channel_name)
    overrides = {
        "donor": args.donor_channel_name,
        "acceptor": args.acceptor_channel_name,
        "laser": args.laser_channel_name,
    }
    for role, override in overrides.items():
        if override and normalize_token(override) in normalized:
            return role
    for role, aliases in CHANNEL_ALIASES.items():
        if any(alias in normalized for alias in aliases):
            return role
    return None


def expand_before_after_roles(files: Sequence[RawNd2File], args: argparse.Namespace) -> Dict[str, RawNd2File]:
    sorted_files = sorted(files, key=lambda item: (item.seq, item.path.name))
    role_files: Dict[str, RawNd2File] = {}
    donor_pending: List[RawNd2File] = []
    acceptor_pending: List[RawNd2File] = []

    for frame_file in sorted_files:
        seq_role = SEQ_ROLE_MAP.get(frame_file.seq)
        if seq_role is not None:
            if seq_role in role_files:
                logging.warning(
                    "Duplicate Seq role %s: keeping %s, ignoring %s",
                    seq_role,
                    role_files[seq_role].path,
                    frame_file.path,
                )
            else:
                role_files[seq_role] = frame_file

    if set(SEQ_ROLE_MAP.values()).issubset(role_files):
        return role_files

    for frame_file in sorted_files:
        if frame_file.role in {"donor_before", "donor_after", "acceptor_before", "acceptor_after", "laser"}:
            if frame_file.role in role_files:
                logging.warning("Duplicate role %s: keeping %s, ignoring %s", frame_file.role, role_files[frame_file.role].path, frame_file.path)
            else:
                role_files[frame_file.role] = frame_file
        elif frame_file.role == "donor":
            donor_pending.append(frame_file)
        elif frame_file.role == "acceptor":
            acceptor_pending.append(frame_file)

    if donor_pending:
        if len(donor_pending) >= 1 and "donor_before" not in role_files:
            role_files["donor_before"] = donor_pending[0]
        if len(donor_pending) >= 2 and "donor_after" not in role_files:
            role_files["donor_after"] = donor_pending[-1]
    if acceptor_pending:
        if len(acceptor_pending) >= 1 and "acceptor_before" not in role_files:
            role_files["acceptor_before"] = acceptor_pending[0]
        if len(acceptor_pending) >= 2 and "acceptor_after" not in role_files:
            role_files["acceptor_after"] = acceptor_pending[-1]

    for role, item in zip(args.role_order, sorted_files):
        role_files.setdefault(role, item)
    return role_files


def discover_raw_measurement_dirs(input_root: Path) -> List[Tuple[str, Path]]:
    subdirs = sorted(p for p in input_root.iterdir() if p.is_dir())
    measurement_dirs = [(p.name, p) for p in subdirs if any(p.rglob("*.nd2"))]
    if any(input_root.glob("*.nd2")):
        measurement_dirs.insert(0, (input_root.name, input_root))
    if not measurement_dirs:
        raise ValueError(f"No ND2 files found under {input_root}.")
    return measurement_dirs


def parse_raw_nd2_files(
    measurement_dir: Path,
    args: argparse.Namespace,
) -> Dict[str, RawNd2File]:
    files: List[RawNd2File] = []
    for path in sorted(measurement_dir.rglob("*.nd2")):
        match = ND2_FILE_RE.match(path.name)
        if not match:
            logging.warning("Skipping ND2 with unexpected name: %s", path)
            continue
        channel = match.group("channel")
        role = infer_channel_role(channel, args)
        if role is None:
            role = "unknown"
        frame = RawNd2File(
            path=path,
            timestamp=match.group("timestamp"),
            token=match.group("token"),
            channel=channel,
            role=role,
            seq=int(match.group("seq")),
        )
        files.append(frame)

    if not files:
        return {}
    return expand_before_after_roles(files, args)


def load_fov_map(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    if path is None:
        return {}
    metadata: Dict[str, Dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            raw_fov = row.get("fov") or row.get("FOV") or row.get("xy") or ""
            if not raw_fov:
                continue
            fov = str(raw_fov).strip()
            if fov.isdigit():
                fov = f"{int(fov):03d}"
            metadata[fov] = {
                "well": row.get("well") or row.get("Well") or "",
                "iptg_label": (
                    row.get("iptg_label")
                    or row.get("IPTG")
                    or row.get("iptg")
                    or row.get("label")
                    or ""
                ),
                "condition": row.get("condition") or row.get("name") or "",
            }
    return metadata


def load_nested_cli_list(value: Optional[str], option_name: str) -> Optional[List[List[object]]]:
    if value is None:
        return None
    candidate = Path(value)
    text = candidate.read_text() if candidate.exists() else value
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{option_name} must be valid JSON or a path to a JSON file.") from exc
    if not isinstance(parsed, list):
        raise ValueError(f"{option_name} must be a list or nested list.")
    if not parsed:
        return []
    if all(not isinstance(item, list) for item in parsed):
        return [list(parsed)]
    if not all(isinstance(item, list) for item in parsed):
        raise ValueError(f"{option_name} must not mix scalar entries and sublists.")
    return [list(item) for item in parsed]


def select_measurement_spec(
    specs: Optional[List[List[object]]],
    measurement_index: int,
    measurement_count: int,
    option_name: str,
) -> Optional[List[object]]:
    if specs is None:
        return None
    if not specs:
        return []
    if len(specs) == measurement_count:
        return specs[measurement_index]
    if len(specs) == 1:
        return specs[0]
    raise ValueError(
        f"{option_name} has {len(specs)} sublists but {measurement_count} measurement "
        "directories were found. Provide one sublist per subdirectory, or one sublist to reuse."
    )


def metadata_for_measurement_position(
    *,
    fov: str,
    position_number: int,
    position_count: int,
    fov_map: Mapping[str, Mapping[str, str]],
    fovs_per_well: int,
    labels: Optional[Sequence[object]],
) -> Tuple[int, str, str, str]:
    if fov in fov_map:
        entry = fov_map[fov]
        well = entry.get("well", "")
        iptg = entry.get("iptg_label", "")
        condition = entry.get("condition", "") or iptg or well or "Multiplexed"
        return 0, well, iptg, condition

    if fovs_per_well <= 0:
        raise ValueError("Raw ND2 mode needs a positive FOVs-per-well value for each measurement.")
    if position_count % fovs_per_well != 0:
        raise ValueError(
            f"{position_count} inferred FOVs is not divisible by {fovs_per_well} FOVs per well."
        )
    group_idx = (position_number - 1) // fovs_per_well
    well_count = position_count // fovs_per_well
    clean_labels = [str(value) for value in labels] if labels else []
    if clean_labels and len(clean_labels) != well_count:
        raise ValueError(
            f"{len(clean_labels)} labels were provided, but {position_count} FOVs at "
            f"{fovs_per_well} FOV/well imply {well_count} wells."
        )

    if clean_labels:
        label = clean_labels[group_idx]
        return group_idx + 1, label, label, label
    well = f"Well{group_idx + 1:02d}"
    iptg = ""
    condition = well
    return group_idx + 1, well, iptg, condition


def read_nd2_position_stack(
    path: Path,
    position_index: int,
    frames_per_fov: Optional[int] = None,
) -> np.ndarray:
    try:
        import nd2
    except ImportError as exc:  # pragma: no cover - depends on local env
        raise ImportError(
            "Raw ND2 mode requires the 'nd2' package. Install/update conda-environment.yaml."
        ) from exc

    with nd2.ND2File(str(path)) as nd_file:
        array = np.asarray(nd_file.asarray())
        axes = list(nd_file.sizes.keys())

    if array.ndim != len(axes):
        axes = axes[-array.ndim :]

    indexers: List[object] = [slice(None)] * array.ndim
    if "P" in axes:
        axis = axes.index("P")
        if position_index < 0 or position_index >= array.shape[axis]:
            raise IndexError(
                f"Position index {position_index} outside ND2 position count {array.shape[axis]} for {path}"
            )
        indexers[axis] = position_index
    elif "T" in axes and frames_per_fov is not None:
        axis = axes.index("T")
        start = position_index * frames_per_fov
        stop = start + frames_per_fov
        if stop > array.shape[axis]:
            raise IndexError(
                f"Position {position_index + 1} needs frames {start}:{stop}, "
                f"but {path} has {array.shape[axis]} time frames."
            )
        indexers[axis] = slice(start, stop)
    elif position_index != 0:
        raise IndexError(
            f"{path} has no ND2 position axis and no frame count for position slicing."
        )

    for singleton_axis in ("C", "Z"):
        if singleton_axis in axes:
            axis = axes.index(singleton_axis)
            indexers[axis] = 0

    array = array[tuple(indexers)]
    kept_axes = [
        axis_name
        for axis_name, selector in zip(axes, indexers)
        if isinstance(selector, slice)
    ]
    array = np.squeeze(array)
    kept_axes = [axis_name for axis_name in kept_axes if axis_name not in {"P", "C", "Z"}]

    if array.ndim == 2:
        array = array[np.newaxis, ...]
    elif "T" in kept_axes:
        time_axis = kept_axes.index("T")
        array = np.moveaxis(array, time_axis, 0)
    elif array.ndim == 3:
        # A single-channel/position ND2 stack commonly arrives as T/Y/X.
        pass
    else:
        while array.ndim > 3:
            array = array[0]
        if array.ndim == 2:
            array = array[np.newaxis, ...]

    if array.ndim != 3:
        raise ValueError(f"Expected a T/Y/X stack from {path}, got shape {array.shape}")
    return array.astype(np.float32, copy=False)


def trim_stack(stack: np.ndarray, count: int, role: str, fov: str) -> np.ndarray:
    if stack.shape[0] < count:
        raise ValueError(f"FOV {fov} {role} has {stack.shape[0]} frames, expected at least {count}.")
    return stack[:count]


def make_cellpose_input(job: RawFovJob, args: argparse.Namespace) -> None:
    donor_file = job.role_files["donor_before"]
    donor_stack = trim_stack(
        read_nd2_position_stack(donor_file.path, job.position_index, args.donor_before_count),
        args.donor_before_count,
        "donor_before",
        job.fov,
    )
    image = np.mean(donor_stack, axis=0).astype(np.float32)
    tifffile.imwrite(job.cellpose_frame_path, image)


def nd2_sizes(path: Path) -> Mapping[str, int]:
    try:
        import nd2
    except ImportError as exc:  # pragma: no cover - depends on local env
        raise ImportError(
            "Raw ND2 mode requires the 'nd2' package. Install/update conda-environment.yaml."
        ) from exc
    with nd2.ND2File(str(path)) as nd_file:
        return {axis: int(size) for axis, size in nd_file.sizes.items()}


def raw_role_frame_counts(args: argparse.Namespace) -> Dict[str, int]:
    return {
        "donor_before": args.donor_before_count,
        "acceptor_before": args.acceptor_before_count,
        "laser": args.laser_count,
        "acceptor_after": args.acceptor_after_count,
        "donor_after": args.donor_after_count,
    }


def infer_position_count(role_files: Mapping[str, RawNd2File], args: argparse.Namespace) -> int:
    expected_counts = raw_role_frame_counts(args)
    inferred: List[int] = []

    for role, raw_file in role_files.items():
        if role not in expected_counts:
            continue
        sizes = nd2_sizes(raw_file.path)
        if "P" in sizes:
            count = sizes["P"]
            if sizes.get("T", expected_counts[role]) < expected_counts[role]:
                raise ValueError(
                    f"{raw_file.path} has only {sizes.get('T', 1)} time frames for {role}; "
                    f"expected {expected_counts[role]} per FOV."
                )
        else:
            total_frames = sizes.get("T", 1)
            per_fov = expected_counts[role]
            if total_frames % per_fov != 0:
                raise ValueError(
                    f"{raw_file.path} has {total_frames} frames for {role}, not divisible by "
                    f"{per_fov} frames/FOV."
                )
            count = total_frames // per_fov
        inferred.append(count)

    if not inferred:
        raise ValueError("Could not infer FOV count from ND2 files.")
    if len(set(inferred)) != 1:
        raise ValueError(f"ND2 files disagree on FOV count: {inferred}")
    return inferred[0]


def fov_label_for_position(position_index: int, args: argparse.Namespace, fov_map: Mapping[str, Mapping[str, str]]) -> str:
    if args.fov_labels and position_index < len(args.fov_labels):
        return str(args.fov_labels[position_index])
    label = f"{position_index + 1:03d}"
    if label in fov_map:
        return label
    return f"pos{position_index + 1:03d}"


def fovs_per_well_for_measurement(
    values: Optional[Sequence[int]],
    measurement_index: int,
    measurement_count: int,
    fallback: int,
) -> int:
    if values:
        if len(values) == measurement_count:
            value = values[measurement_index]
        elif len(values) == 1:
            value = values[0]
        else:
            raise ValueError(
                "--fovs-per-well-by-measurement must contain either one value or one value "
                "per sorted measurement subdirectory."
            )
    else:
        value = fallback
    if value <= 0:
        raise ValueError("FOVs per well must be positive in raw ND2 mode.")
    return int(value)


def build_raw_jobs(args: argparse.Namespace) -> List[RawFovJob]:
    cellpose_dir = ensure_dir(args.output_root / "02_cellpose_raw_input")
    fov_map = load_fov_map(args.fov_map)
    measurement_dirs = discover_raw_measurement_dirs(args.input_root)
    measurement_labels = load_nested_cli_list(args.measurement_labels_json, "--measurement-labels-json")
    jobs: List[RawFovJob] = []

    for measurement_index, (measurement_name, measurement_dir) in enumerate(measurement_dirs):
        role_files = parse_raw_nd2_files(measurement_dir, args)
        if not role_files:
            logging.warning("No usable ND2 files found in %s", measurement_dir)
            continue
        safe_measurement = re.sub(r"[^A-Za-z0-9_.-]+", "_", measurement_name)
        missing_roles = [
            role
            for role in ("donor_before", "donor_after", "acceptor_before", "acceptor_after")
            if role not in role_files
        ]
        if missing_roles:
            raise ValueError(
                f"{measurement_name}: missing required ND2 roles: {', '.join(missing_roles)}. "
                "Use channel names or --role-order to disambiguate the five files."
            )
        if "laser" not in role_files and args.roi_mode == "auto-laser":
            raise ValueError(f"{measurement_name}: --roi-mode auto-laser requires a laser ND2 file.")

        position_count = infer_position_count(role_files, args)
        fovs_per_well = fovs_per_well_for_measurement(
            args.fovs_per_well_by_measurement,
            measurement_index,
            len(measurement_dirs),
            args.fovs_per_well,
        )
        labels_for_measurement = select_measurement_spec(
            measurement_labels,
            measurement_index,
            len(measurement_dirs),
            "--measurement-labels-json",
        )
        if args.fov_labels and len(args.fov_labels) != position_count:
            logging.warning(
                "%s: %d FOV labels provided for %d ND2 positions.",
                measurement_name,
                len(args.fov_labels),
                position_count,
            )

        for position_index in range(position_count):
            fov = fov_label_for_position(position_index, args, fov_map)
            well_index, well, iptg_label, condition = metadata_for_measurement_position(
                fov=fov,
                position_number=position_index + 1,
                position_count=position_count,
                fov_map=fov_map,
                fovs_per_well=fovs_per_well,
                labels=labels_for_measurement,
            )
            frame_path = cellpose_dir / f"{safe_measurement}_{fov}_donor_pre_mean.tif"
            mask_path = frame_path.with_name(frame_path.stem + "_cp_masks.tif")
            job = RawFovJob(
                measurement_name=measurement_name,
                measurement_dir=measurement_dir,
                position_index=position_index,
                fov=fov,
                well_index=well_index,
                well=well,
                iptg_label=iptg_label,
                condition=condition,
                role_files=role_files,
                cellpose_frame_path=frame_path,
                cellpose_mask_path=mask_path,
            )
            make_cellpose_input(job, args)
            jobs.append(job)

    if not jobs:
        raise RuntimeError("No raw ND2 FOV jobs were prepared.")
    return jobs


def read_imagej_roi_mask(roi_path: Path, shape: Tuple[int, int]) -> np.ndarray:
    try:
        from roifile import ImagejRoi
    except ImportError as exc:  # pragma: no cover - depends on local env
        raise ImportError("Reading ImageJ ROI files requires the 'roifile' package.") from exc

    roi = ImagejRoi.fromfile(str(roi_path))
    coords = roi.coordinates()
    if coords is None or len(coords) < 3:
        left = int(getattr(roi, "left", 0))
        right = int(getattr(roi, "right", left + 1))
        top = int(getattr(roi, "top", 0))
        bottom = int(getattr(roi, "bottom", top + 1))
        coords = np.array(
            [[left, top], [right, top], [right, bottom], [left, bottom]],
            dtype=np.float32,
        )
    coords = np.asarray(coords)
    if coords.ndim == 3:
        coords = coords[0]
    rr, cc = draw.polygon(coords[:, 1], coords[:, 0], shape=shape)
    mask = np.zeros(shape, dtype=bool)
    mask[rr, cc] = True
    return mask


RAW_ROI_CACHE: Dict[Tuple[str, Tuple[int, int], str], Tuple[np.ndarray, np.ndarray]] = {}


def prompt_raw_roi_masks(
    job: RawFovJob,
    preview_image: np.ndarray,
    shape: Tuple[int, int],
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    cache_key = (job.measurement_name, shape, args.roi_mode)
    if cache_key in RAW_ROI_CACHE:
        return RAW_ROI_CACHE[cache_key]
    if imagej is None:
        raise ImportError("ROI prompting requires pyimagej. Use --roi-mode auto-laser as a non-GUI fallback.")

    roi_dir = ensure_dir(args.output_root / "rois" / re.sub(r"[^A-Za-z0-9_.-]+", "_", job.measurement_name))
    preview_path = roi_dir / "roi_prompt_preview.tif"
    bleach_roi_path = roi_dir / "bleached.roi"
    unbleached_roi_path = roi_dir / "unbleached.roi"
    tifffile.imwrite(preview_path, preview_image.astype(np.float32, copy=False))

    ask_unbleached = args.roi_mode == "prompt-two"
    ij = init_imagej(args.imagej_distribution)
    macro = f"""
run("Close All");
open("{path_for_macro(preview_path)}");
run("Enhance Contrast", "saturated=0.35");
run("Brightness/Contrast...");
setTool("oval");
waitForUser("Bleached ROI", "Draw the bleached laser area for {job.measurement_name}.\\nThen press OK.");
roiManager("Reset");
roiManager("Add");
roiManager("Select", 0);
roiManager("Save", "{path_for_macro(bleach_roi_path)}");
"""
    if ask_unbleached:
        macro += f"""
setTool("oval");
waitForUser("Unbleached ROI", "Draw the unbleached comparison area for {job.measurement_name}.\\nThen press OK.");
roiManager("Reset");
roiManager("Add");
roiManager("Select", 0);
roiManager("Save", "{path_for_macro(unbleached_roi_path)}");
"""
    macro += """
run("Close All");
"""
    ij.py.run_macro(macro)

    bleach_mask = read_imagej_roi_mask(bleach_roi_path, shape)
    if ask_unbleached:
        ring_mask = read_imagej_roi_mask(unbleached_roi_path, shape)
    else:
        ring_mask = morphology.binary_dilation(
            bleach_mask,
            morphology.disk(max(1, int(args.roi_buffer_px))),
        ) & ~bleach_mask
    RAW_ROI_CACHE[cache_key] = (bleach_mask, ring_mask)
    return bleach_mask, ring_mask


def bleach_mask_from_laser(laser_stack: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    projection = np.max(laser_stack, axis=0)
    if not np.isfinite(projection).any() or np.nanmax(projection) == np.nanmin(projection):
        raise ValueError("Laser stack has no contrast for bleach-mask detection.")

    threshold = np.nanpercentile(projection, args.laser_threshold_percentile)
    candidate = projection >= threshold
    if int(candidate.sum()) < 20:
        threshold = filters.threshold_otsu(projection)
        candidate = projection >= threshold

    candidate = morphology.remove_small_objects(candidate.astype(bool), min_size=20)
    candidate = morphology.binary_closing(candidate, morphology.disk(3))
    labels = measure.label(candidate)
    props = measure.regionprops(labels)
    if not props:
        raise ValueError("Could not derive a bleach mask from the laser channel.")
    largest = max(props, key=lambda prop: prop.area)
    return labels == largest.label


def make_region_masks(
    job: RawFovJob,
    laser_stack: Optional[np.ndarray],
    image_shape: Tuple[int, int],
    preview_image: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    if args.roi_mode in {"prompt-ring", "prompt-two"}:
        return prompt_raw_roi_masks(job, preview_image, image_shape, args)
    if args.laser_roi is not None:
        bleach_mask = read_imagej_roi_mask(args.laser_roi, image_shape)
    elif laser_stack is not None:
        bleach_mask = bleach_mask_from_laser(laser_stack, args)
    else:
        raise ValueError(f"No bleach ROI or laser channel available for {job.measurement_name} {job.fov}.")

    ring_mask = morphology.binary_dilation(
        bleach_mask,
        morphology.disk(max(1, int(args.roi_buffer_px))),
    ) & ~bleach_mask
    return bleach_mask, ring_mask


def background_pixels(
    cell_mask: np.ndarray,
    bleach_mask: np.ndarray,
    args: argparse.Namespace,
) -> np.ndarray:
    occupied = morphology.binary_dilation(
        cell_mask.astype(bool),
        morphology.disk(max(1, int(args.background_exclusion_radius))),
    )
    protected = morphology.binary_dilation(bleach_mask.astype(bool), morphology.disk(2))
    bg_mask = ~(occupied | protected)
    if int(bg_mask.sum()) < 100:
        return ~occupied
    return bg_mask


def scalar_background(stack: np.ndarray, bg_mask: np.ndarray) -> float:
    values: List[float] = []
    for frame in stack:
        pixels = frame[bg_mask]
        if pixels.size >= 50:
            values.append(float(np.nanmedian(pixels)))
        else:
            values.append(float(np.nanpercentile(frame, 10)))
    return float(np.nanmedian(values)) if values else 0.0


def rolling_correct_stack(
    stack: np.ndarray,
    bg_mask: np.ndarray,
    radius: int,
) -> Tuple[np.ndarray, float]:
    corrected = []
    bg_values = []
    for frame in stack:
        background = rolling_ball(frame.astype(np.float32, copy=False), radius=radius)
        corrected.append(np.maximum(frame - background, 0))
        pixels = background[bg_mask]
        bg_values.append(float(np.nanmedian(pixels if pixels.size else background)))
    return np.stack(corrected, axis=0), float(np.nanmedian(bg_values))


def average_stack_region(stack: np.ndarray, region_mask: np.ndarray) -> float:
    if int(region_mask.sum()) == 0:
        return np.nan
    return float(np.nanmean(stack[:, region_mask]))


def classify_cell_area(
    cell_mask: np.ndarray,
    bleach_mask: np.ndarray,
    ring_mask: np.ndarray,
    min_overlap: float,
) -> Tuple[str, float, float]:
    area = max(1, int(cell_mask.sum()))
    bleach_fraction = float(np.logical_and(cell_mask, bleach_mask).sum() / area)
    ring_fraction = float(np.logical_and(cell_mask, ring_mask).sum() / area)
    if bleach_fraction >= min_overlap:
        return "Bleached", bleach_fraction, ring_fraction
    if ring_fraction >= min_overlap:
        return "UnbleachedRing", bleach_fraction, ring_fraction
    return "Outside", bleach_fraction, ring_fraction


def load_raw_quant_stacks(
    job: RawFovJob,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], Dict[str, int]]:
    donor_before = trim_stack(
        read_nd2_position_stack(job.role_files["donor_before"].path, job.position_index, args.donor_before_count),
        args.donor_before_count,
        "donor_before",
        job.fov,
    )
    donor_after = trim_stack(
        read_nd2_position_stack(job.role_files["donor_after"].path, job.position_index, args.donor_after_count),
        args.donor_after_count,
        "donor_after",
        job.fov,
    )
    acceptor_before = trim_stack(
        read_nd2_position_stack(job.role_files["acceptor_before"].path, job.position_index, args.acceptor_before_count),
        args.acceptor_before_count,
        "acceptor_before",
        job.fov,
    )
    acceptor_after = trim_stack(
        read_nd2_position_stack(job.role_files["acceptor_after"].path, job.position_index, args.acceptor_after_count),
        args.acceptor_after_count,
        "acceptor_after",
        job.fov,
    )
    laser_stack = None
    if "laser" in job.role_files:
        laser_stack = trim_stack(
            read_nd2_position_stack(job.role_files["laser"].path, job.position_index, args.laser_count),
            args.laser_count,
            "laser",
            job.fov,
        )
    counts = {
        "donor_before_n": int(donor_before.shape[0]),
        "donor_after_n": int(donor_after.shape[0]),
        "acceptor_before_n": int(acceptor_before.shape[0]),
        "acceptor_after_n": int(acceptor_after.shape[0]),
        "laser_n": 0 if laser_stack is None else int(laser_stack.shape[0]),
    }
    return donor_before, donor_after, acceptor_before, acceptor_after, laser_stack, counts


def quantify_raw_job(job: RawFovJob, args: argparse.Namespace) -> pd.DataFrame:
    mask = load_mask(job.cellpose_mask_path)
    donor_before, donor_after, acc_before, acc_after, laser_stack, counts = load_raw_quant_stacks(job, args)

    image_shape = donor_before.shape[1:]
    for name, stack in {
        "donor_after": donor_after,
        "acceptor_before": acc_before,
        "acceptor_after": acc_after,
    }.items():
        if stack.shape[1:] != image_shape:
            raise ValueError(f"{job.measurement_name} {job.fov}: {name} shape {stack.shape[1:]} != {image_shape}")
    if mask.shape != image_shape:
        raise ValueError(f"{job.measurement_name} {job.fov}: mask shape {mask.shape} != {image_shape}")

    preview_image = np.max(laser_stack, axis=0) if laser_stack is not None else np.mean(donor_before, axis=0)
    bleach_mask, ring_mask = make_region_masks(job, laser_stack, image_shape, preview_image, args)
    cell_union = mask > 0
    bg_mask = background_pixels(cell_union, bleach_mask, args)

    if args.background_mode == "manual":
        if args.background_donor is None or args.background_acceptor is None:
            raise ValueError("--background-mode manual requires --background-donor and --background-acceptor.")
        donor_bg = float(args.background_donor)
        acceptor_bg = float(args.background_acceptor)
        donor_before_corr = donor_before
        donor_after_corr = donor_after
        acc_before_corr = acc_before
        acc_after_corr = acc_after
        subtract_scalar = True
    elif args.background_mode == "rolling-ball":
        donor_before_corr, donor_bg_a = rolling_correct_stack(donor_before, bg_mask, args.rolling_ball_radius)
        donor_after_corr, donor_bg_b = rolling_correct_stack(donor_after, bg_mask, args.rolling_ball_radius)
        acc_before_corr, acceptor_bg_a = rolling_correct_stack(acc_before, bg_mask, args.rolling_ball_radius)
        acc_after_corr, acceptor_bg_b = rolling_correct_stack(acc_after, bg_mask, args.rolling_ball_radius)
        donor_bg = float(np.nanmedian([donor_bg_a, donor_bg_b]))
        acceptor_bg = float(np.nanmedian([acceptor_bg_a, acceptor_bg_b]))
        subtract_scalar = False
    else:
        donor_bg = scalar_background(np.concatenate([donor_before, donor_after], axis=0), bg_mask)
        acceptor_bg = scalar_background(np.concatenate([acc_before, acc_after], axis=0), bg_mask)
        donor_before_corr = donor_before
        donor_after_corr = donor_after
        acc_before_corr = acc_before
        acc_after_corr = acc_after
        subtract_scalar = True

    segmentation_image = np.mean(donor_before, axis=0)
    props = {prop.label: prop for prop in regionprops(mask)}
    records: List[Dict[str, object]] = []

    for label in sorted(value for value in np.unique(mask) if value > 0):
        cell = mask == label
        prop = props.get(int(label))
        if prop is None:
            continue

        pixels = segmentation_image[cell]
        seg_std = float(np.nanstd(pixels))
        quality_ratio = float(np.nanmean(pixels) / seg_std) if seg_std else np.nan
        if np.isnan(quality_ratio) or quality_ratio < args.min_quality_ratio:
            continue

        area_class, bleach_fraction, ring_fraction = classify_cell_area(
            cell,
            bleach_mask,
            ring_mask,
            args.min_area_overlap,
        )

        don_bb = average_stack_region(donor_before, cell)
        don_ab = average_stack_region(donor_after, cell)
        acc_bb = average_stack_region(acc_before, cell)
        acc_ab = average_stack_region(acc_after, cell)

        if subtract_scalar:
            don_bb_bg = don_bb - donor_bg
            don_ab_bg = don_ab - donor_bg
            acc_bb_bg = acc_bb - acceptor_bg
            acc_ab_bg = acc_ab - acceptor_bg
        else:
            don_bb_bg = average_stack_region(donor_before_corr, cell)
            don_ab_bg = average_stack_region(donor_after_corr, cell)
            acc_bb_bg = average_stack_region(acc_before_corr, cell)
            acc_ab_bg = average_stack_region(acc_after_corr, cell)

        fret = np.nan
        if don_ab_bg != 0:
            fret = 100.0 * (don_ab_bg - don_bb_bg) / don_ab_bg

        bleached = np.nan
        if acc_bb_bg != 0:
            bleached = 100.0 * (acc_bb_bg - acc_ab_bg) / acc_bb_bg

        fac = np.nan if np.isnan(bleached) else (100.0 - bleached) / 100.0
        denominator = don_ab_bg - (fac * don_bb_bg if not np.isnan(fac) else 0.0)
        corr_fret = np.nan if denominator == 0 else 100.0 * (don_ab_bg - don_bb_bg) / denominator
        ratio_metric = np.nan if don_ab_bg == 0 else acc_bb_bg / don_ab_bg
        bleached_pass = False if np.isnan(bleached) else bool(bleached >= args.min_bleach_percent)
        analysis_included = (
            area_class == "UnbleachedRing"
            or (area_class == "Bleached" and bleached_pass)
        )

        records.append(
            {
                "run_id": args.run_id,
                "measurement": job.measurement_name,
                "fov": job.fov,
                "position_index": int(job.position_index + 1),
                "well_index": int(job.well_index),
                "well": job.well,
                "iptg_label": job.iptg_label,
                "condition": job.condition,
                "cell_label": int(label),
                "area_class": area_class,
                "analysis_included": bool(analysis_included),
                "in_bleach_roi": area_class == "Bleached",
                "in_unbleached_ring": area_class == "UnbleachedRing",
                "bleach_overlap_fraction": bleach_fraction,
                "ring_overlap_fraction": ring_fraction,
                "area_px": int(cell.sum()),
                "major_axis_px": float(prop.major_axis_length),
                "minor_axis_px": float(prop.minor_axis_length),
                "angle_deg": float(np.degrees(prop.orientation)),
                "quality_ratio": quality_ratio,
                "DonBB": don_bb,
                "DonAB": don_ab,
                "AccBB": acc_bb,
                "AccAB": acc_ab,
                "BG_don": donor_bg,
                "BG_acc": acceptor_bg,
                "BG_pixel_count": int(bg_mask.sum()),
                "DonBB_BG": don_bb_bg,
                "DonAB_BG": don_ab_bg,
                "AccBB_BG": acc_bb_bg,
                "ACCAB_BG": acc_ab_bg,
                "FRET": fret,
                "BleachedPercent": bleached,
                "Fac": fac,
                "CorrFRET": corr_fret,
                "Ratio": ratio_metric,
                "BleachedPass": bleached_pass,
                "background_mode": args.background_mode,
                **counts,
                "cellpose_frame_path": str(job.cellpose_frame_path),
                "cellpose_mask_path": str(job.cellpose_mask_path),
                "measurement_dir": str(job.measurement_dir),
                "donor_before_file": str(job.role_files["donor_before"].path),
                "donor_after_file": str(job.role_files["donor_after"].path),
                "acceptor_before_file": str(job.role_files["acceptor_before"].path),
                "acceptor_after_file": str(job.role_files["acceptor_after"].path),
                "laser_file": str(job.role_files["laser"].path) if "laser" in job.role_files else "",
            }
        )

    return pd.DataFrame(records)


def save_raw_workbook(results_by_measurement: Mapping[str, pd.DataFrame], args: argparse.Namespace) -> Path:
    workbook_path = args.output_xlsx or (args.output_root / "multiplex_results.xlsx")
    ensure_dir(workbook_path.parent)
    used_sheet_names: set[str] = set()
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        wrote_sheet = False
        for measurement_name, df in results_by_measurement.items():
            sheet = sanitize_sheet_name(measurement_name, used_sheet_names)
            if df.empty:
                pd.DataFrame([{"message": "No valid quantified ROIs"}]).to_excel(
                    writer,
                    index=False,
                    sheet_name=sheet,
                )
            else:
                df.to_excel(writer, index=False, sheet_name=sheet)
            wrote_sheet = True
        if not wrote_sheet:
            pd.DataFrame([{"message": "No raw ND2 results generated"}]).to_excel(
                writer,
                index=False,
                sheet_name="NoData",
            )
    logging.info("Saved raw multiplex workbook to %s", workbook_path)
    return workbook_path


def run_raw_nd2_pipeline(args: argparse.Namespace) -> None:
    ensure_dir(args.output_root)
    jobs = build_raw_jobs(args)
    logging.info("Prepared %d raw ND2 FOV jobs.", len(jobs))

    if not args.skip_cellpose:
        run_cellpose_cli(
            cellpose_input=args.output_root / "02_cellpose_raw_input",
            model=args.cellpose_model,
            diameter=args.cellpose_diameter,
            channels=args.cellpose_channels,
        )
    else:
        logging.info("Skipping Cellpose execution.")

    if args.skip_measurement:
        logging.info("Measurement step skipped; raw ND2 pipeline finished after preparation.")
        return

    csv_dir = ensure_dir(args.output_root / "csv")
    results_by_measurement: Dict[str, List[pd.DataFrame]] = {}
    for job in jobs:
        if not job.cellpose_mask_path.exists():
            logging.warning("Missing Cellpose mask for %s FOV %s", job.measurement_name, job.fov)
            continue
        try:
            df = quantify_raw_job(job, args)
        except Exception as exc:
            logging.warning("Skipping %s FOV %s: %s", job.measurement_name, job.fov, exc)
            continue
        if not df.empty:
            results_by_measurement.setdefault(job.measurement_name, []).append(df)

    combined: Dict[str, pd.DataFrame] = {}
    for measurement_name, frames in results_by_measurement.items():
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        combined[measurement_name] = df
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", measurement_name)
        df.to_csv(csv_dir / f"{safe_name}_results.csv", index=False, decimal=args.csv_decimal)

    save_raw_workbook(combined, args)


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
    if imagej is None:
        raise ImportError(
            "pyimagej is required for legacy TIFF/ImageJ mode. "
            "Install the provided conda environment first."
        )
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
    dest_bleached_parent = ensure_dir(config.registered_bleached_root / measurement.name)
    dest_bleached_dir = ensure_dir(dest_bleached_parent / fov)
    dest_unbleached_parent = ensure_dir(config.registered_unbleached_root / measurement.name)
    dest_unbleached_dir = ensure_dir(dest_unbleached_parent / fov)
    roi_mask_path = (dest_dir / "bleach_roi_mask.tif")
    roi_crop_path = (dest_dir / "bleach_roi_crop.roi")
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
// Compute buffered crop rectangle around the bleaching ROI (100 px buffer).
getSelectionBounds(x_roi, y_roi, w_roi, h_roi);
getSelectionCoordinates(xpoints, ypoints);
n_points = xpoints.length;
orig_x = x_roi;
orig_y = y_roi;
crop_x = orig_x - 100;
crop_y = orig_y - 100;
crop_w = w_roi + 200;
crop_h = h_roi + 200;
// Clamp to image bounds
if (crop_x < 0) crop_x = 0;
if (crop_y < 0) crop_y = 0;
img_w = getWidth();
img_h = getHeight();
if (crop_x + crop_w > img_w) crop_w = img_w - crop_x;
if (crop_y + crop_h > img_h) crop_h = img_h - crop_y;
// Shift ROI coordinates into crop coordinate system
for (i = 0; i < n_points; i++) {{
    xpoints[i] = xpoints[i] - crop_x;
    ypoints[i] = ypoints[i] - crop_y;
}}
// Perform the crop using the buffered rectangle
makeRectangle(crop_x, crop_y, crop_w, crop_h);
run("Crop");
run("Enhance Contrast", "saturated=0.35");
resetMinAndMax();
// Recreate bleaching ROI in cropped coordinates
makePolygon(xpoints, ypoints, n_points);
roiManager("Reset");
roiManager("Add");
roiManager("Select", 0);
// Save ROI in crop coordinates and a binary mask image
roi_path = "{path_for_macro(roi_crop_path)}";
roiManager("Save", roi_path);
run("Create Mask");
mask_path = "{path_for_macro(roi_mask_path)}";
saveAs("Tiff", mask_path);
close();
// Save full cropped sequence
    dest = "{path_for_macro(dest_dir)}";
    if (!File.exists(dest)) {{
        File.makeDirectory(dest);
    }}
    print("Saving cropped sequence to " + dest + "/");
    run("Image Sequence... ", "select=[" + dest + "/] dir=[" + dest + "/] format=TIFF name={fov}use");
    print("Saved cropped sequence to " + dest + "/");
    list = getFileList(dest + "/");
    print("Verified " + list.length + " files in destination.");
// Create and save bleached-only stack (inside ROI)
    destBleached = "{path_for_macro(dest_bleached_dir)}";
    if (!File.exists(destBleached)) {{
        File.makeDirectory(destBleached);
    }}
    print("Saving bleached-only sequence to " + destBleached + "/");
    roiManager("Select", 0);
    run("Clear Outside", "stack");
    run("Image Sequence... ", "select=[" + destBleached + "/] dir=[" + destBleached + "/] format=TIFF name={fov}use");
    run("Undo");
// Create and save unbleached-only stack (outside ROI within crop)
    destUnbleached = "{path_for_macro(dest_unbleached_dir)}";
    if (!File.exists(destUnbleached)) {{
        File.makeDirectory(destUnbleached);
    }}
    print("Saving unbleached-only sequence to " + destUnbleached + "/");
    roiManager("Select", 0);
    run("Make Inverse");
    run("Clear Outside", "stack");
    run("Image Sequence... ", "select=[" + destUnbleached + "/] dir=[" + destUnbleached + "/] format=TIFF name={fov}use");
    run("Undo");
    run("Make Inverse");
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
    files = stack_frame_files(stack_dir)
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
    frames = stack_frame_files(stack_dir)
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
    # Load ROI mask (bleached area) saved during registration, if available.
    roi_mask: Optional[np.ndarray] = None
    roi_mask_path = stack_dir / "bleach_roi_mask.tif"
    if roi_mask_path.exists():
        try:
            roi_mask = tifffile.imread(roi_mask_path)
            if roi_mask.ndim == 3:
                roi_mask = roi_mask[0]
            roi_mask = roi_mask.astype(bool)
            if roi_mask.shape != mask.shape:
                logging.warning(
                    "ROI mask shape mismatch for %s-%s: roi %s, mask %s",
                    measurement.name,
                    fov,
                    roi_mask.shape,
                    mask.shape,
                )
                roi_mask = None
        except Exception as exc:
            logging.warning(
                "Failed to load ROI mask for %s-%s from %s: %s",
                measurement.name,
                fov,
                roi_mask_path,
                exc,
            )
            roi_mask = None

    props = {prop.label: prop for prop in regionprops(mask)}
    labels = sorted(label for label in np.unique(mask) if label > 0)
    records: List[Dict[str, object]] = []
    bg_don, bg_acc = config.backgrounds
    
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
        std_series = series.std(axis=1, ddof=0)

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

        # Determine if this cell lies inside the bleaching ROI (In/OUT by geometry).
        in_bleach_roi: Optional[bool] = None
        if roi_mask is not None:
            cy, cx = prop.centroid  # (row, col) = (y, x)
            cy_i = int(round(cy))
            cx_i = int(round(cx))
            if 0 <= cy_i < roi_mask.shape[0] and 0 <= cx_i < roi_mask.shape[1]:
                in_bleach_roi = bool(roi_mask[cy_i, cx_i])
            else:
                in_bleach_roi = False

        def avg_range(start: int, end: int) -> float:
            sl = slice(start - 1, end)
            return float(np.nanmean(mean_series[sl]))

        # Corrected Frame Indices (1-based)
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
            "InBleachROI": in_bleach_roi,
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
        logging.warning("No valid ROIs were quantified in this measurement.")
        return pd.DataFrame()
    df = pd.DataFrame(all_records)
    return df


def save_measurement_results(df: pd.DataFrame, measurement: Measurement, config: PipelineConfig) -> None:
    if df.empty:
        return
    
    ensure_dir(config.output_root)
    
    # Create filename based on measurement name
    csv_path = config.output_root / f"{measurement.name}_results.csv"
    xlsx_path = config.output_root / f"{measurement.name}_results.xlsx"
    
    df.to_csv(csv_path, index=False, decimal=config.csv_decimal)
    logging.info("Saved CSV for %s to %s", measurement.name, csv_path)
    
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Results")
    logging.info("Saved Excel for %s to %s", measurement.name, xlsx_path)


# --------------------------------------------------------------------------------------
# Main orchestration
# --------------------------------------------------------------------------------------


def build_fov_contexts(
    measurements: List[Measurement],
    config: PipelineConfig,
) -> List[FovContext]:
    contexts: List[FovContext] = []
    
    # In Multiplex mode, we don't track global index for IPTG assignment
    # We just process what is given.
    
    num_measurements = len(measurements)
    for m_idx, measurement in enumerate(measurements):
        is_last_measurement = (m_idx == num_measurements - 1)
        # Check if the last measurement contains only the two background FOVs
        is_background_only = is_last_measurement and len(measurement.fovs) == 2
        
        if is_background_only:
            logging.info(
                "Detected that %s contains only 2 FOVs (background FOVs).",
                measurement.name
            )
        
        for idx, fov in enumerate(measurement.fovs, start=1):
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
            
            # Multiplex logic: Mixed cells, but mark background FOVs if detected
            if is_background_only:
                iptg_value = "BLANK"
            else:
                iptg_value = "Multiplexed"
            
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
    args.input_root = args.input_root.resolve()
    args.output_root = args.output_root.resolve()
    if args.output_xlsx is not None:
        args.output_xlsx = args.output_xlsx.resolve()

    if should_use_raw_nd2(args):
        logging.info("Running raw ND2 multiplex mode.")
        run_raw_nd2_pipeline(args)
        return

    if args.background_donor is None or args.background_acceptor is None:
        raise ValueError(
            "Legacy TIFF mode requires --background-donor and --background-acceptor "
            "(or use --mode raw-nd2 with --background-mode auto)."
        )

    # If multiplex flag is on, we enforce multiplex behavior
    # (though this script is designed for it anyway)
    
    config = PipelineConfig(
        input_root=args.input_root.resolve(),
        output_root=args.output_root.resolve(),
        registered_root=(args.output_root / "01_registered").resolve(),
        registered_bleached_root=(args.output_root / "01_registered_bleached").resolve(),
        registered_unbleached_root=(args.output_root / "01_registered_unbleached").resolve(),
        cellpose_input=(args.output_root / "02_cellpose_input").resolve(),
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
        multiplex=args.multiplex or True, # Always True for this script variant
    )

    ensure_dir(config.output_root)
    ensure_dir(config.registered_root)
    ensure_dir(config.registered_bleached_root)
    ensure_dir(config.registered_unbleached_root)
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

    # MULTIPLEX CHANGE: We construct contexts for all, but we need to PROCESS them per measurement.
    # Actually, better to just build and process per measurement loop.
    
    logging.info("Starting Multiplex Analysis (Per-Measurement Output)")
    
    # Cellpose can still be run in batch if we want, OR per measurement.
    # Since Cellpose CLI takes a directory, and we put all frames in 'cellpose_input',
    # it is more efficient to run Cellpose ONCE on the whole folder, then separate the results.
    
    # 1. Extract frames for ALL measurements (to prepare for batch Cellpose)
    all_contexts = build_fov_contexts(measurements, config)
    
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

    # 2. Quantify and Save PER MEASUREMENT
    # We can filter 'all_contexts' by measurement name
    
    for measurement in measurements:
        logging.info("Processing results for %s...", measurement.name)
        
        # Filter contexts for this measurement
        m_contexts = [ctx for ctx in all_contexts if ctx.measurement.name == measurement.name]
        
        if not m_contexts:
            logging.warning("No FOVs found for %s", measurement.name)
            continue
            
        df = process_measurements(m_contexts, config)
        save_measurement_results(df, measurement, config)

    logging.info("Multiplex processing complete.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pragma: no cover
        logging.error("Pipeline failed: %s", exc)
        sys.exit(1)
