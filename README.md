# KD-FRET Pipeline

Unified Python workflow that replicates the legacy KD-FRET processing chain written by Gabriele Malengo
(ImageJ Plugins 1–3 and Excel/Python Plugins 4–9). It:

1. Registers and crops every measurement stack through an ImageJ macro.
2. Extracts Cellpose input frames, optionally runs the Cellpose CLI, and reuses
   the generated ROI masks.
3. Computes all donor/acceptor statistics, FRET metrics, and filters exactly as
   the legacy spreadsheets did.
4. Writes both an appendable CSV (with configurable decimal separator) and a
   ready-to-open Excel workbook containing every ROI record.

---

## Prerequisites

1. **Conda environment**  
   Create it once (Python 3.11) using the provided spec:
   ```bash
   conda env create -f conda-environment.yaml
   conda activate kd-fret
   ```
   The environment bundles `pyimagej`, `nd2`, `tifffile`, `openpyxl`, `cellpose`,
   `roifile`, and the scientific stack.

2. **ImageJ / Fiji**  
   - Default: `pyimagej` downloads a headless Fiji bundle automatically.
   - If you want a specific distribution, provide `--imagej-distribution`.

3. **Input folder layout**  
   - Measurements live under `measurementX` directories (`measurement1`, `measurement2`, ...).
   - Each measurement contains `xyNN` stacks or files that encode `xyNN` in their names.
   - Registered stacks contain 54 frames saved as `xyNNuse0001.tif` … `xyNNuse0054.tif`.

---

## Running the pipeline

### Legacy exported TIFF workflow

```bash
python kd_fret_pipeline.py \
  --input-root /path/to/measurements \
  --output-root /path/to/output \
  --background-donor 170 \
  --background-acceptor 135 \
  --iptg-concentrations 3 10 50 100 200 \
  --fovs-per-iptg 20 \
  --csv-decimal , \
  --results-xlsx /path/to/output/fret_results.xlsx
```

### Raw ND2 multiplex workflow

Use `kd_fret_multiplex.py` when each measurement is a subdirectory containing
the five raw ND2 acquisition stacks. Role assignment is by the trailing `Seq`
digit only; channel substrings in the filename are for logging, not for deciding
which file is donor, acceptor, or laser. mCherry is the acceptor channel. Map:
`Seq0000=acceptor_before`, `Seq0001=donor_before`, `Seq0002=laser`,
`Seq0003=donor_after`, `Seq0004=acceptor_after`. FOVs are inferred from the ND2
position axis, or from the total frame count divided by the expected frames per
FOV when no position axis is present. The default raw frame counts are 2
acceptor-before frames, 20 donor-before frames, 4 laser/bleach frames, 18
donor-after frames, and 2 acceptor-after frames per FOV.

| Seq | Role | Frames/FOV | Analysis use |
| --- | --- | ---: | --- |
| `Seq0000` | `acceptor_before` / mCherry before bleach | 2 | `AccBB`, last 2 frames |
| `Seq0001` | `donor_before` | 20 | Cellpose frame 5; `DonBB`, last 3 frames |
| `Seq0002` | laser / acceptor bleaching | 4 | stack alignment continuity and ROI context |
| `Seq0003` | `donor_after` | 18 | `DonAB`, first 3 frames |
| `Seq0004` | `acceptor_after` / mCherry after bleach | 2 | ROI prompt image; `AccAB`, first 2 frames |

With `--nd2-alignment sift` (default), the script writes a flat virtual TIFF
multiplex under `<output-root>/nd2_ij_source/`, runs Fiji linear stack alignment
(SIFT), prompts for one global laser ROI on an aligned post-bleach mCherry image,
then writes cropped aligned stacks under `01_registered/`. Cellpose receives the
fifth donor pre-bleach frame from the cropped aligned stack, matching the legacy
TIFF workflow. Pass `--laser-roi` to skip the prompt and reuse a saved ROI. Use
`--nd2-alignment none` only as a non-SIFT fallback; in that mode `--roi-mode`
controls prompt-ring, prompt-two, or auto-laser ROI handling.

To restart from saved aligned/cropped stacks, keep the same input/output roots
and add `--skip-registration`. The raw SIFT restart path reuses
`<output-root>/01_registered/<measurement>/xyNN/`. Add `--skip-cellpose` as well
when the matching masks already exist under
`<output-root>/03_cellpose_raw_output/<measurement>/`.

```bash
python kd_fret_multiplex.py \
  --mode raw-nd2 \
  --input-root /path/to/day \
  --output-root /path/to/output \
  --fovs-per-well-by-measurement 5 10 10 5 \
  --background-mode auto
```

Raw mode writes one workbook, `<output-root>/multiplex_results.xlsx`, with one
sheet per measurement subdirectory. `--fovs-per-well-by-measurement` is one
integer per sorted subdirectory; the pipeline infers the well count from
`total FOVs / FOVs per well`. If labels are omitted, wells are named `Well01`,
`Well02`, etc. Add `--measurement-labels-json` only when you want custom well or
condition labels. The output rows include `BG_don`, `BG_acc`, `BG_pixel_count`,
and `background_mode`.

Workflow summary:

1. The script discovers all measurements/FOVs and initializes ImageJ unless
   `--skip-registration` is reusing saved registered stacks.
2. If no ROI is supplied, it aligns one stack, shows the post-bleach mCherry
   frame in ImageJ, and waits for you to draw the laser ROI (oval tool selected
   by default, brightness/contrast dialog open).
3. Registers and crops every stack, pauses 10 seconds between measurements (to
   mimic Plugin 1), prepares Cellpose frames, and optionally runs Cellpose.
4. Loads each registered stack plus the mirrored `03_cellpose_raw_output/` mask,
   applies the ROI quality filter, computes Don/Acc metrics, FRET, CorrFRET,
   ratio, etc.
5. Writes one `.xlsx` workbook with one sheet per measurement subdirectory and
   per-measurement CSV files under `<output-root>/csv/`.

---

## Multiplex CLI Reference

### Paths and outputs
| Flag | Description |
| --- | --- |
| `--input-root PATH` | Required. Folder containing measurement subdirectories or raw ND2 files. |
| `--output-root PATH` | Root for intermediates and final outputs. |
| `--mode auto\|raw-nd2\|legacy-tiff` | Auto-detect raw ND2 input or force a mode. |
| `--output-xlsx PATH` | Raw ND2 workbook path, default `<output-root>/multiplex_results.xlsx`. |

### Registration & ROI
| Flag | Description |
| --- | --- |
| `--laser-roi PATH` | Reuse an existing ImageJ ROI file instead of prompting. |
| `--nd2-alignment sift\|none` | Raw ND2 alignment mode; `sift` is the deployment/default path. |
| `--roi-mode prompt-ring\|prompt-two\|auto-laser` | Only used with `--nd2-alignment none`. |
| `--skip-registration` | Reuse existing registered stacks; raw SIFT expects `01_registered/<measurement>/xyNN/`. |
| `--imagej-distribution STR` | Pass custom ImageJ/Fiji distribution string to `imagej.init()`. |
| `--sequence-start INT` | Legacy TIFF start index passed to `File.openSequence` (default 5). |
| `--nd2-align-frame-start INT` | Raw ND2 SIFT start index; keep the default `1`. |

### Cellpose
| Flag | Description |
| --- | --- |
| `--cellpose-frame-index INT` | Frame extracted per FOV for Cellpose (default 5). |
| `--cellpose-model STR` | Cellpose pretrained model (default `bact_fluor_cp3`). |
| `--cellpose-diameter FLOAT` | Expected object diameter (default 25). |
| `--cellpose-channels CH1 CH2` | Channel indices for Cellpose (default `0 0`). |
| `--cellpose-review-outputs` | Raw ND2 only: also save Cellpose PNG/outlines/ROI review files; slower on old CPUs. |
| `--skip-cellpose` | Skip launching Cellpose; assumes `_cp_masks.tif` already exist. |

### Measurement logic
| Flag | Description |
| --- | --- |
| `--background-mode manual\|auto\|rolling-ball` | Raw ND2 background correction mode. |
| `--background-donor FLOAT` | Donor background value, required for manual background and legacy TIFF mode. |
| `--background-acceptor FLOAT` | Acceptor background value, required for manual background and legacy TIFF mode. |
| `--fovs-per-well-by-measurement ...` | FOVs per well for each sorted measurement subdirectory, or one value reused for all. |
| `--measurement-labels-json JSON_OR_PATH` | Optional nested labels per sorted measurement subdirectory. |
| `--fov-map PATH` | Optional CSV mapping FOVs to well/condition metadata. |
| `--min-quality-ratio FLOAT` | Mean/StdDev threshold for ROI inclusion (default 0.8). |
| `--min-bleach-percent FLOAT` | Bleach filter for bleached cells, default 70.0. |
| `--run-id STR` | Run label stored in outputs (default: timestamp). |
| `--csv-decimal "."\|","` | Decimal separator for CSV write (default `.`). |
| `--skip-measurement` | Stops after preparing registration/cellpose inputs. |

### Misc
| Flag | Description |
| --- | --- |
| `--log-level LEVEL` | `DEBUG`, `INFO`, `WARNING`, or `ERROR` (default `INFO`). |

---

## Outputs

```
<output-root>/
├── 01_registered/MeXX/xyYY/*.tif        # Cropped, aligned stacks
├── 02_cellpose_raw_input/MeXX/*.tif     # Raw ND2 Cellpose inputs, mirrored by measurement
├── 03_cellpose_raw_output/MeXX/*_cp_masks.tif
├── csv/*_results.csv                    # Per-measurement CSV files
└── multiplex_results.xlsx               # One sheet per measurement subdirectory
```

Every row in the CSV/XLSX corresponds to a single Cellpose ROI and contains:

- Measurement metadata (name, FOV, well, condition, run ID)
- ROI geometry (area, axes, angle), quality ratio, background-corrected signals
- FRET, Bleached%, Fac, CorrFRET, ratio (AccBB_BG / DonAB_BG)
- References to the registered stack directory and mask path

---

## Troubleshooting tips

- **ROI prompt image looks pixelated** – use the Brightness/Contrast dialog shown by the macro. No resizing occurs; if it still looks blocky, inspect the corresponding `01_registered/` stack to confirm acquisition quality.
- **IPTG values appear as 30 instead of 3** – open the CSV with the same decimal separator you used via `--csv-decimal`. Set `,` for European Excel.
- **Negative or huge FRET values** – verify donor/acceptor background numbers match the acquisition, and inspect the stored `DonBB/DonAB` values in the CSV/XLSX.
- **Missing Cellpose masks** – ensure you didn’t run with `--skip-cellpose` before the `_cp_masks.tif` files exist.

For any other issue, re-run with `--log-level DEBUG` and check the console plus the ImageJ window output; the macros now print detailed status messages around registration, saving, and ROI handling.
