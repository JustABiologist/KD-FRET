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
   Create it once (Python 3.10) using the provided spec:
   ```bash
   conda env create -f conda-environment.yaml
   conda activate kd-fret
   ```
   The environment already bundles `pyimagej`, `tifffile`, `openpyxl`, `cellpose`,
   and the scientific stack.

2. **ImageJ / Fiji**  
   - Default: `pyimagej` downloads a headless Fiji bundle automatically.
   - If you want a specific distribution, provide `--imagej-distribution`.

3. **Input folder layout**  
   - Measurements live under `measurementX` directories (`measurement1`, `measurement2`, ...).
   - Each measurement contains `xyNN` stacks or files that encode `xyNN` in their names.
   - Registered stacks contain 54 frames saved as `xyNNuse0001.tif` … `xyNNuse0054.tif`.

---

## Running the pipeline

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

Workflow summary:

1. The script discovers all measurements/FOVs, initializes ImageJ (unless reused
   ROI and `--skip-registration`).
2. If no ROI is supplied, it aligns one stack, shows the post-bleach mCherry
   frame in ImageJ, and waits for you to draw the laser ROI (oval tool selected
   by default, brightness/contrast dialog open).
3. Registers+ crops every stack, pauses 10 seconds between measurements (to mimic
   Plugin 1), prepares Cellpose frames, and optionally runs Cellpose.
4. Loads each registered stack plus the `_cp_masks.tif` mask, applies the ROI
   quality filter, computes Don/Acc metrics, FRET, CorrFRET, ratio, etc.
5. Appends results to the CSV (respecting the requested decimal separator) and
   writes the same table to an `.xlsx` workbook.

---

## CLI reference

### Paths and outputs
| Flag | Description |
| --- | --- |
| `--input-root PATH` | **Required.** Folder containing `measurementX` subfolders. |
| `--output-root PATH` | Root for all intermediate outputs (`01_registered`, `02_cellpose_input`, results). |
| `--existing-results PATH` | Append to an existing CSV instead of starting fresh. |
| `--results-xlsx PATH` | Optional Excel path (default `<output-root>/fret_results.xlsx`). |

### Registration & ROI
| Flag | Description |
| --- | --- |
| `--laser-roi PATH` | Reuse an existing ImageJ ROI file instead of prompting. |
| `--skip-registration` | Assume registered stacks already exist (requires `--laser-roi`). |
| `--imagej-distribution STR` | Pass custom ImageJ/Fiji distribution string to `imagej.init()`. |
| `--sequence-start INT` | Frame index passed to `File.openSequence` (default 5). |

### Cellpose
| Flag | Description |
| --- | --- |
| `--cellpose-frame-index INT` | Frame extracted per FOV for Cellpose (default 5). |
| `--cellpose-model STR` | Cellpose pretrained model (default `bact_fluor_cp3`). |
| `--cellpose-diameter FLOAT` | Expected object diameter (default 25). |
| `--cellpose-channels CH1 CH2` | Channel indices for Cellpose (default `0 0`). |
| `--skip-cellpose` | Skip launching Cellpose; assumes `_cp_masks.tif` already exist. |

### Measurement logic
| Flag | Description |
| --- | --- |
| `--background-donor FLOAT` | Donor background value used for subtraction. |
| `--background-acceptor FLOAT` | Acceptor background value. |
| `--min-quality-ratio FLOAT` | Mean/StdDev threshold for ROI inclusion (default 0.8). |
| `--min-bleach-percent FLOAT` | Bleach filter (default 69.6, matching DaySummaryIN). |
| `--iptg-concentrations ...` | List of concentrations consumed sequentially across the full dataset (last measurement’s final two FOVs automatically become `BLANK`). |
| `--fovs-per-iptg INT` | FOV count per concentration block (default 20). |
| `--default-group STR` | Fallback `IN/OUT` label if measurement name doesn’t specify one. |
| `--run-id STR` | Run label stored in outputs (default: timestamp). |
| `--csv-decimal "."|","` | Decimal separator for CSV read/write (default `.`). |
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
├── 02_cellpose_input/*.tif/PNG/ZIP      # Cellpose inputs + masks/ROIs
├── fret_results.csv                     # Appendable CSV (decimal separator configurable)
└── fret_results.xlsx                    # Identical data in Excel format (sheet “Results”)
```

Every row in the CSV/XLSX corresponds to a single Cellpose ROI and contains:

- Measurement metadata (name, IN/OUT, FOV, IPTG concentration, run ID)
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