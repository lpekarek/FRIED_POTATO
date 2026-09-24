"""
Paramater_optimization_analysis.py

Compares batch FRIED POTATO analyses against a manual (TOMATO) ground-truth
analysis and ranks the parameter sets.

Expected inputs
---------------
MANUAL_FILE : CSV exported from the TOMATO table ('Save results table').
              Columns used: filename, step number, mean force [pN],
              extension step start [nm], Step length [nm].

BATCH files : total_results_*.csv files written by the folder analysis.
              Columns used: filename, orientation, step number/step #,
              Fc, step start, step end, step length.
              Malformed rows (known exporter bug: duplicated/partial step
              assignments) are repaired by truncation/padding.

Optionally, a parameters_<timestamp>.txt next to each total_results_<timestamp>.csv
(the file written by export_settings) is parsed, so the ranking table shows
the actual parameter values behind each score.
"""

import csv
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

# ======================================================================
# CONFIGURATION - adjust everything here
# ======================================================================

FOLDER = r"D:\OT_data_raw\03_OT_data_Groningen\test_parameters\Total_results_files"
MANUAL_FILE = os.path.join(FOLDER, "Manual.csv")
BATCH_PATTERN = os.path.join(FOLDER, "total_results_*.csv")
OUTPUT_FILE = os.path.join(FOLDER, "parameter_ranking.csv")

# --- Matching / scoring tolerances ---
FORCE_TOL_PN = 0.5       # manual vs batch mean-force agreement window
DIST_TOL_NM = 1.0         # manual vs batch step-length agreement window
MATCH_TOL_NM = 5.0        # max |difference in step start| for two steps
                          # to be considered the SAME physical event

# --- Scoring weights (per curve, score starts at 1000, higher = better) ---
STEP_COUNT_WEIGHT = 100   # penalty per missing OR extra step
PARAM_PENALTY_FACTOR = 10 # multiplier for parameter deviations beyond tolerance
BASE_SCORE = 1000

# --- Diagnostics ---
NEAR_DUPLICATE_NM = 2.0  # two steps in the same curve whose step start is
                          # closer than this = suspected exporter duplicate

# Tokens exported by the fitting code that mean "no usable value"
BAD_TOKENS = {'none', 'nan', '#n/a', '', 'null', 'na'}

# Processing suffixes found in TOMATO/manual filenames, e.g.
# "..._fw_curve1_smooth_20260922-123934_downsamplesd"
SUFFIX_RE = re.compile(r"(_smooth_\d{8}-\d{6})?(downsample[sd]+)?$", flags=re.IGNORECASE)


# ======================================================================
# GENERIC HELPERS
# ======================================================================

def clean_filename(fn):
    """Normalize curve IDs so manual and batch files can be matched."""
    if fn is None or (isinstance(fn, float) and np.isnan(fn)):
        return fn
    fn = str(fn).strip()
    fn = SUFFIX_RE.sub("", fn)
    return fn


def _clean_token(v):
    """Map failed-fit placeholders ('None', 'nan', '#N/A') to real NaN."""
    if isinstance(v, str):
        s = v.strip()
        return np.nan if s.lower() in BAD_TOKENS else s
    return v


def robust_read_csv(filepath):
    """
    Read a TOMATO/FRIED-POTATO export that may contain corrupted rows:
      - rows with MORE fields than the header (duplicate/partial write
        concatenation bug) -> truncated, first N fields kept
      - rows with FEWER fields -> padded with NaN
      - failure tokens -> converted to NaN
    Returns (DataFrame, n_anomalies).
    """
    with open(filepath, 'r', encoding='utf-8-sig', newline='') as fh:
        rows = list(csv.reader(fh))

    if not rows:
        return pd.DataFrame(), 0

    header = [h.strip().lower() for h in rows[0]]
    ncols = len(header)

    clean_rows, anomalies = [], 0
    for i, row in enumerate(rows[1:], start=2):
        if not row or all(v.strip() == '' for v in row):
            continue                      # skip blank lines
        if len(row) != ncols:
            anomalies += 1
            row = row[:ncols] if len(row) > ncols else row + [''] * (ncols - len(row))
        clean_rows.append([_clean_token(v) for v in row])

    if anomalies:
        print(f"  Note: {os.path.basename(filepath)}: {anomalies} malformed row(s) "
              f"(header has {ncols} fields) - repaired by truncation/padding.")

    return pd.DataFrame(clean_rows, columns=header), anomalies


# ======================================================================
# DATA LOADING
# ======================================================================

def load_steps_from_csv(filepath, source):
    """
    Extract the minimal comparable data from manual OR batch CSVs.

    source: 'manual' or 'batch'
    Returns DataFrame with columns:
        filename, step_number, force (pN), step_start (nm), step_length (nm)
    """
    df, anomalies = robust_read_csv(filepath)
    if df.empty:
        return pd.DataFrame(), anomalies

    if 'filename' not in df.columns:
        raise KeyError(f"'filename' column not found in {filepath}. "
                       f"Available columns: {list(df.columns)}")

    df['filename'] = df['filename'].apply(clean_filename)

    # step-number column is called 'step number' in total_results,
    # 'step #' in the per-file/common-steps exports
    step_col = 'step number' if 'step number' in df.columns else 'step #'
    if step_col not in df.columns:
        raise KeyError(f"No step number column in {filepath}. "
                       f"Available: {list(df.columns)}")

    df['step_number'] = pd.to_numeric(df[step_col], errors='coerce')
    df = df[df['step_number'] > 0].copy()       # step 0 = baseline row

    if source == 'batch':
        force_col, start_col, length_col = 'fc', 'step start', 'step length'
    else:
        force_col, start_col, length_col = 'mean force [pn]', 'extension step start [nm]', 'step length [nm]'

    for col, out in ((force_col, 'force'), (start_col, 'step_start'), (length_col, 'step_length')):
        if col not in df.columns:
            df[out] = np.nan
        else:
            df[out] = pd.to_numeric(df[col], errors='coerce')

    keep = ['filename', 'step_number', 'force', 'step_start', 'step_length']
    if 'orientation' in df.columns and source == 'batch':
        keep.append('orientation')

    df = df[keep]

    # Drop rows where NONE of the comparables exist (pure baseline junk)
    df = df.dropna(subset=['force', 'step_start', 'step_length'], how='all')

    # De-duplicate exact repeats of the same (curve, step number)
    df = df.drop_duplicates(subset=['filename', 'step_number'])

    # Warn about near-duplicate step starts (exporter bug signature)
    for fn, grp in df.groupby('filename'):
        starts = grp['step_start'].dropna().sort_values().values
        if len(starts) > 1:
            gaps = np.diff(starts)
            if (gaps < NEAR_DUPLICATE_NM).any():
                print(f"  Warning: {fn}: {int((gaps < NEAR_DUPLICATE_NM).sum())} "
                      f"near-duplicate step start(s) < {NEAR_DUPLICATE_NM} nm apart "
                      f"- suspected duplicated-step export.")

    return df.sort_values(['filename', 'step_number']).reset_index(drop=True), anomalies


def load_parameters_txt(batch_csv_path):
    """
    Find and parse the parameters_<timestamp>.txt belonging to a
    total_results_<timestamp>.csv (written by export_settings).

    Returns a flat dict of Data processing + Input format parameters.
    """
    m = re.search(r"total_results_(\d{8}-\d{6})\.csv$", os.path.basename(batch_csv_path))
    candidates = []
    if m:
        candidates.append(os.path.join(os.path.dirname(batch_csv_path),
                                       f"parameters_{m.group(1)}.txt"))
    candidates.extend(sorted(glob.glob(os.path.join(
        os.path.dirname(batch_csv_path), "parameters_*.txt"))))

    for cand in candidates:
        if os.path.isfile(cand):
            try:
                return parse_parameters_txt(cand)
            except Exception as e:
                print(f"  Warning: could not parse {cand}: {e}")
                return {}
    return {}


def parse_parameters_txt(txt_path):
    """Parse the three JSON sections written by export_settings()."""
    with open(txt_path, 'r') as fh:
        lines = fh.readlines()

    def section(start_marker, end_markers):
        try:
            start = lines.index(start_marker)
        except ValueError:
            return {}
        end = len(lines)
        for em in end_markers:
            try:
                end = min(end, lines.index(em, start + 1))
            except ValueError:
                pass
        try:
            return json.loads(''.join(lines[start + 1:end]).strip())
        except json.JSONDecodeError:
            return {}

    data_processing = section("Data processing:\n", ["Fitting parameters:\n", "Input format:\n"])
    fitting = section("Fitting parameters:\n", ["Input format:\n"])
    input_format = section("Input format:\n", [])

    params = {f"P_{k}": v for k, v in data_processing.items()}
    params.update({f"F_{k}": v for k, v in fitting.items()})
    params.update({f"I_{k}": v for k, v in input_format.items()})
    return params


# ======================================================================
# SCORING
# ======================================================================

def match_steps(m_steps, b_steps, tol_nm):
    """
    Greedy nearest-neighbour matching of manual and batch steps within one
    curve, based on step start position. Steps closer than tol_nm are
    considered the same physical event.

    Returns (matches, missed, extra) where matches is a list of
    (manual_idx, batch_idx).
    """
    pairs = []
    for mi, ms in m_steps.iterrows():
        if pd.isna(ms['step_start']):
            continue
        for bi, bs in b_steps.iterrows():
            if pd.isna(bs['step_start']):
                continue
            d = abs(ms['step_start'] - bs['step_start'])
            if d <= tol_nm:
                pairs.append((d, mi, bi))

    pairs.sort(key=lambda t: t[0])          # closest pairs first
    used_m, used_b, matches = set(), set(), []
    for d, mi, bi in pairs:
        if mi not in used_m and bi not in used_b:
            used_m.add(mi)
            used_b.add(bi)
            matches.append((mi, bi))

    missed = [i for i in m_steps.index if i not in used_m]
    extra = [i for i in b_steps.index if i not in used_b]
    return matches, missed, extra


def score_curve(m_steps, b_steps,
                force_tol=FORCE_TOL_PN, dist_tol=DIST_TOL_NM,
                match_tol=MATCH_TOL_NM,
                step_weight=STEP_COUNT_WEIGHT,
                param_factor=PARAM_PENALTY_FACTOR):
    """
    Score one curve. Returns dict with score and detailed diagnostics.
    """
    matches, missed, extra = match_steps(m_steps, b_steps, match_tol)

    force_devs, dist_devs = [], []
    for mi, bi in matches:
        mf, bf = m_steps.loc[mi, 'force'], b_steps.loc[bi, 'force']
        ml, bl = m_steps.loc[mi, 'step_length'], b_steps.loc[bi, 'step_length']
        if pd.notna(mf) and pd.notna(bf):
            force_devs.append(abs(mf - bf))
        if pd.notna(ml) and pd.notna(bl):
            dist_devs.append(abs(ml - bl))

    n_missed, n_extra = len(missed), len(extra)

    # Parameter penalty: only deviations BEYOND tolerance count
    excess = 0.0
    for d in force_devs:
        excess += max(0.0, d - force_tol)
    for d in dist_devs:
        excess += max(0.0, d - dist_tol)
    avg_excess = excess / len(matches) if matches else 0.0

    param_penalty = avg_excess * param_factor
    step_penalty = (n_missed + n_extra) * step_weight

    score = max(0.0, BASE_SCORE - step_penalty - param_penalty)

    return {
        'score': round(score, 1),
        'n_manual': len(m_steps),
        'n_batch': len(b_steps),
        'n_matched': len(matches),
        'n_missed': n_missed,
        'n_extra': n_extra,
        'mean_force_dev_pN': round(np.mean(force_devs), 3) if force_devs else np.nan,
        'mean_len_dev_nm': round(np.mean(dist_devs), 3) if dist_devs else np.nan,
    }


def evaluate_batch(manual_df, batch_df, only_orientation=None):
    """
    Evaluate one batch run against the manual analysis over all shared
    curves. Returns the aggregate score and per-metric averages.
    """
    per_curve = {}
    curves = [c for c in manual_df['filename'].unique() if c in set(batch_df['filename'])]

    for curve in curves:
        m = manual_df[manual_df['filename'] == curve]
        b = batch_df[batch_df['filename'] == curve]
        if only_orientation and 'orientation' in b.columns:
            b = b[b['orientation'] == only_orientation]
        if b.empty:
            per_curve[curve] = {'score': 0, 'n_manual': len(m), 'n_batch': 0,
                                'n_matched': 0, 'n_missed': len(m), 'n_extra': 0,
                                'mean_force_dev_pN': np.nan, 'mean_len_dev_nm': np.nan}
            continue
        per_curve[curve] = score_curve(m, b)

    if not per_curve:
        return None, {}

    agg = pd.DataFrame(per_curve).T
    overall = {
        'score': agg['score'].mean(),
        'curves_compared': len(agg),
        'steps_matched': int(agg['n_matched'].sum()),
        'steps_missed': int(agg['n_missed'].sum()),
        'steps_extra': int(agg['n_extra'].sum()),
        'precision': (agg['n_matched'].sum() /
                      max(1, (agg['n_matched'] + agg['n_extra']).sum())),
        'recall': (agg['n_matched'].sum() /
                   max(1, (agg['n_matched'] + agg['n_missed']).sum())),
        'mean_force_dev_pN': agg['mean_force_dev_pN'].astype(float).mean(),
        'mean_len_dev_nm': agg['mean_len_dev_nm'].astype(float).mean(),
    }
    return overall, per_curve


# ======================================================================
# MAIN
# ======================================================================

def main():
    print("=" * 70)
    print("FRIED POTATO -- Parameter Optimization Analysis")
    print("=" * 70)

    # --- Manual (ground truth) ---
    print(f"\nLoading manual analysis: {MANUAL_FILE}")
    try:
        manual_df, _ = load_steps_from_csv(MANUAL_FILE, source='manual')
    except Exception as e:
        sys.exit(f"Failed to load manual file: {type(e).__name__}: {e}")
    if manual_df.empty:
        sys.exit("Manual file contained no usable steps.")

    n_curves = manual_df['filename'].nunique()
    print(f"  -> {len(manual_df)} manual steps across {n_curves} curves.")

    # --- Batch files ---
    files = sorted(glob.glob(BATCH_PATTERN))
    files = [f for f in files if not os.path.basename(f).startswith('Manual')]
    if not files:
        sys.exit(f"No batch files found matching: {BATCH_PATTERN}")
    print(f"\nFound {len(files)} batch result files.")

    results = []
    per_curve_dump = {}

    for f in files:
        name = os.path.basename(f)
        print(f"\nProcessing {name} ...")
        try:
            batch_df, anomalies = load_steps_from_csv(f, source='batch')
        except Exception as e:
            print(f"  SKIP - {type(e).__name__}: {e}")
            continue
        if batch_df.empty:
            print("  SKIP - no usable steps.")
            continue

        overall, per_curve = evaluate_batch(manual_df, batch_df)
        if overall is None:
            print("  SKIP - no curve overlap with manual analysis "
                  "(check filename normalization!).")
            continue

        params = load_parameters_txt(f)
        params['file'] = name
        params['malformed_rows'] = anomalies
        params.update({k: (round(v, 3) if isinstance(v, float) else v)
                       for k, v in overall.items()})
        results.append(params)
        per_curve_dump[name] = per_curve

    if not results:
        sys.exit("No batch file could be scored.")

    results_df = pd.DataFrame(results)

    # Put score + diagnostics first, parameter columns after
    front = ['file', 'score', 'curves_compared', 'steps_matched',
             'steps_missed', 'steps_extra', 'precision', 'recall',
             'mean_force_dev_pN', 'mean_len_dev_nm', 'malformed_rows']
    cols = [c for c in front if c in results_df.columns] + \
           [c for c in results_df.columns if c not in front]
    results_df = results_df[cols].sort_values('score', ascending=False).reset_index(drop=True)

    # Rank column
    results_df.insert(0, 'rank', results_df.index + 1)

    print("\n" + "=" * 70)
    print("RANKING (best to worst)")
    print("=" * 70)
    with pd.option_context('display.max_columns', None, 'display.width', 250):
        print(results_df.head(15).to_string(index=False))

    results_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nFull ranking saved to: {OUTPUT_FILE}")

    # Per-curve detail for the winner
    best = results_df.iloc[0]['file']
    print(f"\nBest parameter set: '{best}'")
    if best in per_curve_dump:
        det = pd.DataFrame(per_curve_dump[best]).T
        print("Per-curve breakdown for the winning run:")
        with pd.option_context('display.max_rows', None, 'display.width', 200):
            print(det.to_string())

    worst_curves = det[det['score'] < BASE_SCORE]
    if len(worst_curves) > 0:
        print(f"\nNote: {len(worst_curves)} curve(s) were scored below "
              f"{BASE_SCORE} even for the best run - consider reviewing "
              f"them individually (possibly noisy or unusual curves).")


if __name__ == "__main__":
    main()