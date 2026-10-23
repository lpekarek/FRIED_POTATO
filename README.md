# FRIED POTATO
**Force-Ramp Improved EDition of Practical Optical Tweezers Analysis TOol**

**Version:** 3.1 (2026-10-23)
**Developers:** Lukáš Pekárek & Stefan Buck & Simon Doll & Leo König 
**Copyright:** © 2026 Lukáš Pekárek & Stefan Buck & Simon Doll

---

## 📖 Overview

FRIED POTATO is an automated Python-based GUI tool for processing Force-Distance
(FD) data from optical tweezers experiments. It detects unfolding events, fits
polymer models (WLC/FJC) to the curve sections, and calculates the mechanical
work applied to the studied structure (e.g., RNA).

![FRIED POTATO](https://github.com/lpekarek/FRIED_POTATO/assets/59534208/05407fc0-f169-4a62-b842-3ba7ccd2fcfe)

It extends the original [POTATO](https://github.com/REMI-HIRI/POTATO) tool with
improved exports, additional analysis modes, and robustness enhancements.

For the full method description of the original tool, refer to the POTATO
manuscript: doi: [10.1101/2021.11.11.468103](https://doi.org/10.1101/2021.11.11.468103)

---

## ✨ What's New in FRIED POTATO (vs. original POTATO)

1. **Extended export columns** — change in contour length per step (`delta Lc`),
   total contour length (`total Lc`), total work (`total W`) and step counts are
   written to the cumulative `total_results_<timestamp>.csv` file.
2. **Reverse Fitting mode** — suitable for long RNA analysis when full unfolding
   has been reached; fits the fully unfolded region first and propagates the
   ds-parameters into the ss-region fits (works in both batch and TOMATO modes).
3. **Savitzky-Golay filtering** — optional alternative to the Butterworth filter;
   labels and parameters switch automatically between polynomial order /
   window length and filter degree / cut-off frequency.
4. **Trap Position-based curve splitting (MultiH5)** — when `Trap position/1X`
   data is present in an h5 file, forward/reverse curves are split using trap
   movement rather than bead distance, avoiding misclassification caused by
   unfolding steps. Falls back to distance automatically for older files.
5. **SVG export** — all figures (batch analysis plots, constant-force plots and
   fitted histograms) can additionally be exported as vector graphics.
6. **Robustness improvements** — automated error logging
   (`error_log_<timestamp>.txt` in the analysis folder), empty/short-array
   validation before orientation checks and derivatives, per-file crash recovery
   (failed fits are recorded as `NaN` rows with a `fit_status` message instead of
   aborting the run), and explicit memory cleanup between files.
7. **Descriptive analysis folder names** — output folders are named with a short
   parameter code, e.g. `Analysis_M0_DER_DS10_BW_FD3_CO1e+04_FM5_ZF2_ZD2_...`,
   encoding input format, step-detection mode, downsampling, filter settings,
   force threshold and z-scores used.
8. **TOMATO total-results export** — the manual TOMATO analysis now writes a
   batch-compatible cumulative `total_results_<timestamp>.csv` plus a per-curve
   `<filename>_fit_data_<timestamp>.csv` containing the model curves.
9. **Parameter load/save** — parameter sets saved as text files (JSON sections)
   can be reloaded via *File → Load Parameters*, restoring all preprocessing,
   fitting and input-format settings.

---

## 🛠️ Dependencies

**Python 3** with the following packages:
`csv`, `glob`, `h5py`, `json`, `lumicks.pylake`, `matplotlib`, `multiprocessing`,
`numpy`, `os`, `pandas`, `pathlib`, `PIL`, `scipy`, `statistics`, `time`, `tkinter`

A conda environment file (`environment.yml`) with pinned package versions is
provided. A standalone Windows executable is also available to avoid any
environment setup.

**Launch (source version):**
```bash
python FRIED_POTATO_GUI.py
```

---

## 📂 Input Formats

| Format | Source | Frequency |
|---|---|---|
| **High Frequency (Piezo Distance)** | `.h5` (LUMICKS C-Trap, piezo tracking) | Read from the h5 file |
| **Low Frequency** | `.h5` (LUMICKS C-Trap) | Calculated from start/stop timestamps in the file |
| **CSV (F/D)** | Two columns, no headers: Force (pN), Distance (µm or nm) | Entered by the user |

Notes:
- CSV data can be forward (unfolding) or reverse (folding) curves. For the script
  to distinguish them, reverse curves must **start with the highest distance**.
- Force channel (Trap 1x / 2x) and distance unit (µm / nm) are selectable in the GUI.
- **ATTENTION:** Switching input format resets all parameters to that format's defaults.

---

## 🖥️ GUI Tabs

### 1. Folder Analysis
Automated batch analysis of all files in a selected folder.

- Select the input format, adjust parameters, then press
  **"Select Folder to Analyse!"** (or *File → Analyse folder*).
- Each curve's orientation is detected automatically; reverse curves are flipped
  internally and labelled `orientation = reverse` in the outputs.
- With **MultiH5** enabled, multi-pull h5 files are split into individual
  forward/reverse curves (named `..._fw_curveN` / `..._rv_curveN`), optionally
  using Trap Position data.
- Optional **Min step length** mode replaces the classic z-score/derivative
  step detection with contour-length-based step filtering.
- Optional **Data Augmentation** duplicates/jitters data to stabilise fits.

**Outputs** (in `Analysis[_RF]_<parameter-code>_<timestamp>/`):
- `*_smooth_*.csv` — down-sampled, filtered FD data
- PNG (+ optional SVG) plots of the FD curve, derivatives and marked steps
- `*_results_*.csv` — steps found per file (force derivative, distance
  derivative, and common steps)
- Fitting plots and parameter tables per curve
- `total_results_<timestamp>.csv` — all steps from all files, combined with the
  fit results (`model_type`, `log_likelihood`, ds/ss contour lengths,
  persistence lengths, stiffnesses, offsets and their standard errors, work in
  pN·nm and k<sub>B</sub>T, `delta Lc`, `total Lc`, `total W`, `fit_status`)
- `parameters_<timestamp>.txt` — all settings used (re-loadable)
- The **Refresh** button displays the latest generated plot and prints
  progress messages in the output window.

### 2. Show Single File
Displays a single FD curve (h5-HF, h5-LF or CSV) twice: once raw/unfiltered and
once down-sampled and filtered. Useful for verifying preprocessing settings and
diagnosing data problems before running a batch analysis. *File → Show h5 file
structure* prints the group/dataset layout of an h5 file to help locate data.

![Tab2_new](https://user-images.githubusercontent.com/77443646/173436271-b313fdc2-908f-47fb-bae0-527c56b25a2f.PNG)

### 3. Constant Force Analysis
Analysis of constant-force markers.

  
![Tab3_new](https://user-images.githubusercontent.com/77443646/173436319-698ef81b-68c2-4957-8cc7-579660421fb6.PNG)

1. **Display Constant Force Data** — time-distance plot with an adjacent
   distance histogram.
2. Enter the number of expected Gaussians plus initial guesses for mean, SD and
   amplitude (comma-separated, lowest to highest).
3. **Fit Constant Force Data** — the histogram is fitted with a sum of Gaussians;
   per-peak areas and population fractions are computed.

**Outputs** (in `Analysis_constantF_<file>_<timestamp>/`): the time-distance
figure, the fitted histogram figure (PNG **and** SVG), the smoothed data, the
histogram distribution, and a parameter table (mean, sigma, amplitude, area,
fraction), each as CSV.

### 4. Advanced Settings
All adjustable parameters, grouped by analysis stage:

  
![Tab4_new](https://user-images.githubusercontent.com/77443646/173436359-56f1a257-432f-4d47-826e-6b4efc288e10.PNG)

**Preprocessing**
- Downsampling rate
- Filter type (Butterworth or Savitzky-Golay), filter degree / polynomial order
- Cut-off frequency / SG window length
- Force threshold (pN)

**Derivative**
- Step *d* — interval used for numerical derivatives
- Data frequency (Hz)

**Statistics**
- Z-scores (force / distance)
- Moving median window size
- SD difference threshold (iteration stopping criterion)

**Fitting**
- Model choice: `WLC+WLC` or `WLC+FJC`
- ds/ss persistence lengths (ds/ss *L*<sub>p</sub>), contour lengths (*L*<sub>c</sub>),
  stretch moduli (*K*<sub>0</sub>) — each with upper/lower bounds
- Force and distance offsets, with bounds

**Export**
Toggles for: Processed FD data · Plot (PNG/SVG) · Steps found · Total results ·
Fitting. Unticking *Fitting* skips all fitting steps and greatly accelerates the
run.

### 5. TOMATO — Manual Analysis
Interactive step-picking for difficult curves, using the same preprocessing and
fitting settings as the batch analysis.

![Tab5_new](https://user-images.githubusercontent.com/77443646/173436395-dca85422-a263-4aed-b87b-f18f7eec3aaf.PNG)

- **Choose folder**, then mark each (un)folding step in order from lowest to
  highest: `s` (or *Set start*) + click · `e` (or *Set end*) + click ·
  `Ctrl+S` (or *Save step*). The last saved step can be removed with
  *Delete step*.
- **Analyze curve** (`Ctrl+F`) performs fitting and work calculations on the
  marked regions — the GUI is unresponsive during this, please be patient.
  Bad fits can be removed from the results table (select row + `Del`).
- Navigate between curves with **←/→** or `a`/`d` (step table resets each time).
- **Save results table** exports the full results table as CSV.
- Per curve, the fitted model curves are exported to
  `<filename>_fit_data_<timestamp>.csv`, and each analysed curve is appended to a
  batch-compatible `total_results_<timestamp>.csv` with the same columns as the
  folder analysis.
- Reverse Fitting is also available in TOMATO via the same tick-box.

---

## ⚠️ Troubleshooting

- **"No file of the selected data type in the folder!"** — check that the folder
  contains `.h5` or `.csv` files matching the selected input format (HF/LF/CSV).
- **Curves below the force threshold are skipped** — check that the correct trap
  (1x/2x) and units (µm/nm) are selected.
- **Crashes / failed fits** — consult `error_log_<timestamp>.txt` inside the
  analysis folder; failed individual fits are logged in the total-results CSV
  with a `fit_status` message instead of stopping the analysis.

---

## 📚 References

- Original POTATO manuscript:
  doi:[10.1101/2021.11.11.468103](https://doi.org/10.1101/2021.11.11.468103)
- LUMICKS C-Trap: https://www.lumicks.com

---

FRIED POTATO — *may your unfoldings be crisp.*
