import pandas as pd
import glob
import os
import re
import numpy as np

# Regex to strip processing suffixes like "_smooth_20260922-123934_downsamplesd"
SUFFIX_RE = re.compile(r"_smooth_\d{8}-\d{6}_downsample\w*", flags=re.IGNORECASE)

def clean_filename(fn):
    """Normalize curve IDs so manual and batch files can be matched."""
    if pd.isna(fn):
        return fn
    fn = str(fn).strip()
    fn = SUFFIX_RE.sub("", fn)
    return fn

def normalize_columns(df):
    """Lowercase and trim column names, fixing 'Filename' vs 'filename'."""
    df.columns = [c.strip().lower() for c in df.columns]
    return df

def load_manual_data(filepath):
    """Loads and cleans the manual analysis CSV."""
    df = pd.read_csv(filepath)
    df = normalize_columns(df)

    if 'filename' not in df.columns:
        raise KeyError(f"'filename' column not found in {filepath}. "
                       f"Columns are: {list(df.columns)}")

    df['filename'] = df['filename'].apply(clean_filename)

    # Make sure 'step number' is numeric before filtering
    df['step number'] = pd.to_numeric(df['step number'], errors='coerce')
    df_steps = df[df['step number'] > 0].copy()

    # Drop duplicated (curve, step) rows (some baseline rows repeat in exports)
    df_steps = df_steps.drop_duplicates(subset=['filename', 'step number'])

    force_cols = ['force step start [pn]', 'force step end [pn]', 'mean force [pn]']
    ext_cols = ['extension step start [nm]', 'extension step end [nm]', 'step length [nm]']

    for col in force_cols + ext_cols:
        if col in df_steps.columns:
            df_steps[col] = pd.to_numeric(df_steps[col], errors='coerce')

    return df_steps

def load_batch_file(filepath):
    """Loads a single batch analysis CSV, cleaned and filtered."""
    df = pd.read_csv(filepath)
    df = normalize_columns(df)

    if 'filename' not in df.columns:
        print(f"Warning: 'filename' column not found in {os.path.basename(filepath)}, skipping.")
        return None

    df['filename'] = df['filename'].apply(clean_filename)
    df['step number'] = pd.to_numeric(df['step number'], errors='coerce')
    df_steps = df[df['step number'] > 0].copy()
    df_steps = df_steps.drop_duplicates(subset=['filename', 'step number'])

    # Ensure comparison columns are numeric
    for col in ['fc', 'step length']:
        if col in df_steps.columns:
            df_steps[col] = pd.to_numeric(df_steps[col], errors='coerce')

    return df_steps

def calculate_score(manual_df, batch_df,
                    force_threshold=2.0,
                    dist_threshold=1.0,
                    step_count_weight=100,
                    param_penalty_factor=10.0):
    """Similarity score between manual and batch analysis. Higher is better."""
    scores = {}
    curves = manual_df['filename'].unique()

    for curve in curves:
        m_rows = manual_df[manual_df['filename'] == curve].sort_values('step number')
        b_rows = batch_df[batch_df['filename'] == curve].sort_values('step number')

        if len(b_rows) == 0:
            scores[curve] = 0
            continue

        # 1. Step Count Penalty
        diff_steps = abs(len(m_rows) - len(b_rows))
        step_penalty = diff_steps * step_count_weight

        # 2. Parameter Matching (index-based pairing of ordered steps)
        param_diff_sum = 0
        matched_count = 0
        min_len = min(len(m_rows), len(b_rows))

        for i in range(min_len):
            m_row = m_rows.iloc[i]
            b_row = b_rows.iloc[i]

            m_force = m_row.get('mean force [pn]', np.nan)
            b_force = b_row.get('fc', np.nan)

            m_dist = m_row.get('step length [nm]', np.nan)
            b_dist = b_row.get('step length', np.nan)

            if pd.notna(m_force) and pd.notna(b_force):
                excess_force = max(0, abs(m_force - b_force) - force_threshold)
                param_diff_sum += excess_force

            if pd.notna(m_dist) and pd.notna(b_dist):
                excess_dist = max(0, abs(m_dist - b_dist) - dist_threshold)
                param_diff_sum += excess_dist

            matched_count += 1

        if matched_count > 0:
            avg_param_penalty = (param_diff_sum / matched_count) * param_penalty_factor
        else:
            avg_param_penalty = 1000

        curve_score = 1000 - step_penalty - avg_param_penalty
        scores[curve] = max(0, curve_score)

    return sum(scores.values()) / len(scores) if scores else 0

def main():
    folder_path = r"D:\OT_data_raw\03_OT_data_Groningen\test_parameters\Total_results_files"
    MANUAL_FILE = os.path.join(folder_path, "Manual.csv")
    BATCH_PATTERN = os.path.join(folder_path, "total_results_*.csv")

    FORCE_TOL_PN = 0.3
    DIST_TOL_NM = 1.0

    print("--- FRIED POTATO Optimization Script ---")
    print(f"Loading manual data from: {MANUAL_FILE}")

    try:
        manual_df = load_manual_data(MANUAL_FILE)
        n_curves = manual_df['filename'].nunique()
        print(f"Loaded {len(manual_df)} steps across {n_curves} curves from manual analysis.")

        files = glob.glob(BATCH_PATTERN)
        if not files:
            print(f"No batch files found matching: {BATCH_PATTERN}")
            return
        print(f"\nFound {len(files)} batch result files.")

        results = []
        for f in files:
            batch_df = load_batch_file(f)
            if batch_df is None or batch_df.empty:
                continue

            # Sanity check: how many manual curves were found in this batch file
            overlap = manual_df['filename'].isin(batch_df['filename']).sum()
            score = calculate_score(
                manual_df, batch_df,
                force_threshold=FORCE_TOL_PN,
                dist_threshold=DIST_TOL_NM
            )

            results.append({
                "file": os.path.basename(f),
                "score": round(score, 2),
                "matched_steps_in_manual": overlap
            })

        results_df = pd.DataFrame(results).sort_values(by='score', ascending=False)

        print("\n--- Ranking Results (Best to Worst) ---")
        print(results_df.to_string(index=False))

        if not results_df.empty:
            print(f"\nRecommended Parameters: Based on '{results_df.iloc[0]['file']}'")

        # Optional: save rankings to CSV
        out_path = os.path.join(folder_path, "parameter_ranking.csv")
        results_df.to_csv(out_path, index=False)
        print(f"\nRankings saved to: {out_path}")

    except Exception as e:
        print(f"Error during execution: {type(e).__name__}: {e}")

if __name__ == "__main__":
    main()