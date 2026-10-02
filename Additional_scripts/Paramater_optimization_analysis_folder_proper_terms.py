import pandas as pd
import glob
import os
import re
import numpy as np
from pathlib import Path
import json

# Strip processing suffixes from filenames
SUFFIX_RE = re.compile(r"_smooth_\d{8}-\d{6}_downsample\w*", flags=re.IGNORECASE)

def clean_filename(fn):
    """Normalize curve IDs so manual and batch files can be matched."""
    if pd.isna(fn):
        return fn
    fn = str(fn).strip()
    fn = SUFFIX_RE.sub("", fn)
    return fn

def normalize_columns(df):
    """Lowercase and trim column names."""
    df.columns = [c.strip().lower() for c in df.columns]
    return df

def load_manual_data(filepath):
    """Load and clean manual analysis CSV."""
    df = pd.read_csv(filepath)
    df = normalize_columns(df)

    if 'filename' not in df.columns:
        raise KeyError(f"'filename' column not found in {filepath}")

    df['filename'] = df['filename'].apply(clean_filename)
    df['step number'] = pd.to_numeric(df['step number'], errors='coerce')

    # Focus on actual unfolding events (step > 0)
    df_steps = df[df['step number'] > 0].copy()
    df_steps = df_steps.drop_duplicates(subset=['filename', 'step number'])

    # Ensure numeric columns
    force_cols = ['mean force [pn]', 'force step start [pn]', 'force step end [pn]']
    ext_cols = ['step length [nm]', 'extension step start [nm]', 'extension step end [nm]']

    for col in force_cols + ext_cols:
        if col in df_steps.columns:
            df_steps[col] = pd.to_numeric(df_steps[col], errors='coerce')

    return df_steps

def load_batch_file(filepath):
    """Load a single batch analysis CSV."""
    df = pd.read_csv(filepath)
    df = normalize_columns(df)

    if 'filename' not in df.columns:
        print(f"Warning: 'filename' column not found in {os.path.basename(filepath)}")
        return None

    df['filename'] = df['filename'].apply(clean_filename)
    df['step number'] = pd.to_numeric(df['step number'], errors='coerce')
    df_steps = df[df['step number'] > 0].copy()
    df_steps = df_steps.drop_duplicates(subset=['filename', 'step number'])

    for col in ['fc', 'step length']:
        if col in df_steps.columns:
            df_steps[col] = pd.to_numeric(df_steps[col], errors='coerce')

    return df_steps

def parse_parameters_txt(filepath):
    """Parse parameters.txt from FRIED POTATO analysis folder."""
    params = {}
    if not os.path.exists(filepath):
        return params

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        lines = content.split('\n')

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if '=' in line:
                parts = line.split('=', 1)
                if len(parts) == 2:
                    params[parts[0].strip()] = parts[1].strip()
            elif ':' in line and not line.startswith('/'):
                parts = line.split(':', 1)
                if len(parts) == 2:
                    params[parts[0].strip()] = parts[1].strip()

    except Exception as e:
        print(f"Warning: Could not parse {filepath}: {e}")

    return params

def find_manual_file(search_root):
    """Search for Manual.csv going up the directory tree."""
    search_paths = [Path(search_root)]
    current = Path(search_root).parent
    for _ in range(5):
        search_paths.append(current)
        current = current.parent

    for path in search_paths:
        candidates = list(path.glob("Manual.csv"))
        if candidates:
            return candidates[0]
    return None

def find_matching_steps(manual_df, batch_df,
                        force_tol=0.3,
                        dist_tol=1.0,
                        method='extension_overlap'):
    """
    Match manual steps to batch steps using greedy assignment.

    Returns:
    --------
    Dictionary with matched pairs and unmatched items per curve.
    """
    curve_matches = {}
    curves = manual_df['filename'].unique()

    for curve in curves:
        m_rows = manual_df[manual_df['filename'] == curve].sort_values('step number').reset_index(drop=True)
        b_rows = batch_df[batch_df['filename'] == curve].sort_values('step number').reset_index(drop=True)

        if len(b_rows) == 0:
            # --- FIXED: include n_manual and n_batch in ALL return paths ---
            curve_matches[curve] = {
                'matched_pairs': [],
                'unmatched_manual': list(range(len(m_rows))),
                'unmatched_batch': [],
                'n_manual': len(m_rows),
                'n_batch': 0
            }
            continue

        # Build extension ranges
        manual_ranges = []
        batch_ranges = []

        for i, row in m_rows.iterrows():
            step_len = row.get('step length [nm]', np.nan)
            manual_ranges.append((i, step_len if pd.notna(step_len) else 0))

        for i, row in b_rows.iterrows():
            step_len = row.get('step length', np.nan)
            batch_ranges.append((i, step_len if pd.notna(step_len) else 0))

        # Greedy matching based on extension center proximity
        matched_m = set()
        matched_b = set()
        matches = []

        # Calculate all pairwise distances
        distances = []
        for mi, (m_idx, m_len) in enumerate(manual_ranges):
            for bi, (b_idx, b_len) in enumerate(batch_ranges):
                if m_len > 0 and b_len > 0:
                    center_diff = abs(m_len - b_len)
                    norm_diff = center_diff / max(m_len, b_len)
                    distances.append((mi, bi, center_diff, norm_diff))

        # Sort by distance (closest first)
        distances.sort(key=lambda x: x[3])

        # Greedy assignment
        for mi, bi, center_diff, norm_diff in distances:
            if mi not in matched_m and bi not in matched_b:
                if center_diff <= dist_tol or norm_diff <= 0.3:  # 30% relative tolerance
                    matches.append((mi, bi, center_diff, norm_diff))
                    matched_m.add(mi)
                    matched_b.add(bi)

        unmatched_manual = [i for i in range(len(m_rows)) if i not in matched_m]
        unmatched_batch = [i for i in range(len(b_rows)) if i not in matched_b]

        curve_matches[curve] = {
            'matched_pairs': matches,
            'unmatched_manual': unmatched_manual,
            'unmatched_batch': unmatched_batch,
            'n_manual': len(m_rows),
            'n_batch': len(b_rows)
        }

    return curve_matches

def calculate_classification_metrics(curve_matches,
                                     force_tol=0.3,
                                     dist_tol=1.0):
    """
    Calculate standard classification metrics for step detection.

    Definitions:
    ------------
    TP (True Positive): Batch step matches manual step within tolerance
    FN (False Negative): Manual step NOT detected by batch
    FP (False Positive): Batch step that doesn't match any manual step
    """
    tp_total = 0
    fp_total = 0
    fn_total = 0

    per_curve_metrics = {}

    for curve, match_info in curve_matches.items():
        # --- FIXED: defensive access so malformed entries can never crash ---
        matched_pairs = match_info.get('matched_pairs', [])
        n_manual = match_info.get('n_manual', len(match_info.get('unmatched_manual', [])))
        n_batch = match_info.get('n_batch', len(match_info.get('unmatched_batch', [])))

        # TP: number of matched pairs
        tp = len(matched_pairs)

        # FN: manual steps not matched
        fn = n_manual - tp

        # FP: batch steps not matched
        fp = n_batch - tp

        tp_total += tp
        fp_total += fp
        fn_total += fn

        # Per-curve precision/recall
        if (tp + fp) > 0:
            precision = tp / (tp + fp)
        else:
            precision = 0.0

        if (tp + fn) > 0:
            recall = tp / (tp + fn)
        else:
            recall = 0.0

        if (precision + recall) > 0:
            f1 = 2 * precision * recall / (precision + recall)
        else:
            f1 = 0.0

        per_curve_metrics[curve] = {
            'tp': tp,
            'fp': fp,
            'fn': fn,
            'precision': precision,
            'recall': recall,
            'f1': f1
        }

    # Aggregate metrics
    if (tp_total + fp_total) > 0:
        precision = tp_total / (tp_total + fp_total)
    else:
        precision = 0.0

    if (tp_total + fn_total) > 0:
        recall = tp_total / (tp_total + fn_total)
    else:
        recall = 0.0

    if (precision + recall) > 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = 0.0

    # Accuracy approximation: proportion of correctly classified "events"
    total_events = tp_total + fp_total + fn_total
    if total_events > 0:
        accuracy_approx = tp_total / total_events
    else:
        accuracy_approx = 1.0

    # MCC: TN is not defined for step detection, use F1 as proxy
    mcc = f1

    return {
        'tp': tp_total,
        'fp': fp_total,
        'fn': fn_total,
        'tn': np.nan,  # Not applicable
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'accuracy': accuracy_approx,
        'specificity': np.nan,  # Not applicable
        'mcc': mcc,
        'per_curve': per_curve_metrics
    }

def calculate_detailed_metrics(manual_df, batch_df,
                               force_tol=0.3,
                               dist_tol=1.0):
    """Calculate all metrics including force/distance deviations."""

    # Get classification metrics
    curve_matches = find_matching_steps(manual_df, batch_df,
                                        force_tol=force_tol,
                                        dist_tol=dist_tol)
    class_metrics = calculate_classification_metrics(curve_matches,
                                                     force_tol=force_tol,
                                                     dist_tol=dist_tol)

    # Add force/distance deviation stats for matched pairs
    force_devs = []
    dist_devs = []

    for curve, match_info in curve_matches.items():
        m_rows = manual_df[manual_df['filename'] == curve].sort_values('step number').reset_index(drop=True)
        b_rows = batch_df[batch_df['filename'] == curve].sort_values('step number').reset_index(drop=True)

        for mi, bi, _, _ in match_info['matched_pairs']:
            m_row = m_rows.iloc[mi]
            b_row = b_rows.iloc[bi]

            m_force = m_row.get('mean force [pn]', np.nan)
            b_force = b_row.get('fc', np.nan)

            m_dist = m_row.get('step length [nm]', np.nan)
            b_dist = b_row.get('step length', np.nan)

            if pd.notna(m_force) and pd.notna(b_force):
                force_devs.append(abs(m_force - b_force))

            if pd.notna(m_dist) and pd.notna(b_dist):
                dist_devs.append(abs(m_dist - b_dist))

    class_metrics['mean_force_dev_pN'] = np.nanmean(force_devs) if force_devs else np.nan
    class_metrics['std_force_dev_pN'] = np.nanstd(force_devs) if force_devs else np.nan
    class_metrics['mean_dist_dev_nm'] = np.nanmean(dist_devs) if dist_devs else np.nan
    class_metrics['std_dist_dev_nm'] = np.nanstd(dist_devs) if dist_devs else np.nan

    return class_metrics

def compute_composite_score(class_metrics,
                            f1_weight=0.4,
                            recall_weight=0.3,
                            precision_weight=0.2,
                            param_weight=0.1):
    """
    Compute composite quality score.
    Higher is better (0-100 scale).
    """
    f1 = class_metrics.get('f1_score', 0)
    recall = class_metrics.get('recall', 0)
    precision = class_metrics.get('precision', 0)

    # Penalize parameter deviations
    mean_force_dev = class_metrics.get('mean_force_dev_pN', np.nan)
    mean_dist_dev = class_metrics.get('mean_dist_dev_nm', np.nan)

    force_penalty = min(1.0, mean_force_dev / 1.0) if not np.isnan(mean_force_dev) else 0
    dist_penalty = min(1.0, mean_dist_dev / 2.0) if not np.isnan(mean_dist_dev) else 0
    param_quality = 1.0 - (force_penalty + dist_penalty) / 2.0

    score = (f1_weight * f1 +
             recall_weight * recall +
             precision_weight * precision +
             param_weight * param_quality)

    return score * 100.0  # Scale to 0-100

def find_analysis_folders(root_path):
    """Find all Analysis_* timestamp folders recursively."""
    root = Path(root_path)
    analysis_folders = []

    for folder in root.rglob("Analysis_*"):
        if folder.is_dir():
            csv_files = list(folder.glob("total_results*.csv"))
            if csv_files:
                analysis_folders.append(folder)

    return sorted(analysis_folders)

def main():
    # Configuration
    ROOT_PATH = r"D:\OT_data_raw\03_OT_data_Groningen\test_parameters"

    # Matching tolerances (tune based on your data quality)
    FORCE_TOL_PN = 0.3      # +/- pN tolerance for "matching" step
    DIST_TOL_NM = 1.0       # +/- nm tolerance for "matching" step

    # Composite score weights (tune based on priorities)
    F1_WEIGHT = 0.4         # F1 score (most important - balance)
    RECALL_WEIGHT = 0.3     # Recall (don't miss steps)
    PRECISION_WEIGHT = 0.2  # Precision (avoid spurious steps)
    PARAM_WEIGHT = 0.1      # Parameter accuracy (secondary)

    print("=" * 80)
    print("FRIED POTATO OPTIMIZER -- Classification Metrics & Step Detection Quality")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  Root path: {ROOT_PATH}")
    print(f"  Force tolerance: ±{FORCE_TOL_PN} pN")
    print(f"  Distance tolerance: ±{DIST_TOL_NM} nm")
    print(f"  Score weights: F1={F1_WEIGHT}, Recall={RECALL_WEIGHT}, Precision={PRECISION_WEIGHT}, Param={PARAM_WEIGHT}")
    print()

    # 1. Load manual analysis
    manual_file = find_manual_file(ROOT_PATH)
    if manual_file is None:
        print(f"Error: Manual.csv not found near {ROOT_PATH}")
        return

    print(f"Manual analysis file: {manual_file}")

    try:
        manual_df = load_manual_data(str(manual_file))
        n_manual_curves = manual_df['filename'].nunique()
        n_manual_steps = len(manual_df)
        print(f"Loaded {n_manual_steps} steps across {n_manual_curves} curves\n")
    except Exception as e:
        print(f"Error loading manual data: {e}")
        return

    # 2. Find all analysis folders
    analysis_folders = find_analysis_folders(ROOT_PATH)
    if not analysis_folders:
        print(f"No Analysis_* folders found under {ROOT_PATH}")
        return

    print(f"Found {len(analysis_folders)} analysis folders:\n")
    for af in analysis_folders:
        print(f"  - {af.name}")
    print()

    # 3. Score each run with classification metrics
    results = []

    for folder in analysis_folders:
        print(f"Processing: {folder.name}")

        # Load parameters
        params_file = folder / "parameters.txt"
        params = parse_parameters_txt(params_file)

        # Load batch results
        csv_files = list(folder.glob("total_results*.csv"))
        if not csv_files:
            print(f"  Warning: No total_results*.csv found")
            continue

        all_batches = []
        for csv_file in csv_files:
            batch_df = load_batch_file(str(csv_file))
            if batch_df is not None and not batch_df.empty:
                all_batches.append(batch_df)

        if not all_batches:
            print(f"  Warning: No valid step data")
            continue

        batch_df = pd.concat(all_batches, ignore_index=True)

        # Calculate comprehensive metrics
        metrics = calculate_detailed_metrics(manual_df, batch_df,
                                            force_tol=FORCE_TOL_PN,
                                            dist_tol=DIST_TOL_NM)

        # Compute composite score
        score = compute_composite_score(metrics,
                                        f1_weight=F1_WEIGHT,
                                        recall_weight=RECALL_WEIGHT,
                                        precision_weight=PRECISION_WEIGHT,
                                        param_weight=PARAM_WEIGHT)

        # Extract key parameters for display
        display_params = {}
        for key in ['threshold', 'smoothing', 'step_size', 'min_step_length',
                    'F_min', 'filter_type', 'downsample_value']:
            if key in params:
                display_params[key] = params[key]

        results.append({
            'analysis_folder': folder.name,
            'folder_path': str(folder),
            'score': round(score, 1),
            'tp': metrics['tp'],
            'fp': metrics['fp'],
            'fn': metrics['fn'],
            'precision': round(metrics['precision'], 3),
            'recall': round(metrics['recall'], 3),
            'f1_score': round(metrics['f1_score'], 3),
            'accuracy': round(metrics['accuracy'], 3),
            'mean_force_dev_pN': round(metrics['mean_force_dev_pN'], 2) if not np.isnan(metrics['mean_force_dev_pN']) else 'N/A',
            'mean_dist_dev_nm': round(metrics['mean_dist_dev_nm'], 2) if not np.isnan(metrics['mean_dist_dev_nm']) else 'N/A',
            'params_json': str(display_params)
        })

        print(f"  → Score: {score:.1f}/100 | F1: {metrics['f1_score']:.3f} | "
              f"P: {metrics['precision']:.3f} | R: {metrics['recall']:.3f} | "
              f"TP:{metrics['tp']} FP:{metrics['fp']} FN:{metrics['fn']}")

    # 4. Rank and display results
    if not results:
        print("\nNo valid results to rank.")
        return

    results_df = pd.DataFrame(results).sort_values(by='score', ascending=False)

    print("\n" + "=" * 80)
    print("RANKING BY CLASSIFICATION PERFORMANCE")
    print("=" * 80)

    def format_params(param_str):
        try:
            d = json.loads(param_str.replace("'", '"'))
            return ", ".join([f"{k}={v}" for k, v in list(d.items())[:4]])
        except:
            return param_str[:50] + "..." if len(param_str) > 50 else param_str

    summary_df = results_df.copy()
    summary_df['params_summary'] = summary_df['params_json'].apply(format_params)

    cols_to_show = ['analysis_folder', 'score', 'f1_score', 'precision', 'recall',
                    'tp', 'fp', 'fn', 'params_summary']
    print(summary_df[cols_to_show].to_string(index=False))

    # 5. Save results
    output_file = Path(ROOT_PATH) / "parameter_optimization_ranking.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nDetailed rankings saved to: {output_file}")

    # Human-readable report
    report_file = Path(ROOT_PATH) / "parameter_optimization_report.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("FRIED POTATO Parameter Optimization Report\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Manual Analysis: {manual_file}\n")
        f.write(f"Total Curves: {n_manual_curves}\n")
        f.write(f"Total Manual Steps: {n_manual_steps}\n")
        f.write(f"Force Tolerance: ±{FORCE_TOL_PN} pN\n")
        f.write(f"Distance Tolerance: ±{DIST_TOL_NM} nm\n\n")
        f.write(f"Score Weights: F1={F1_WEIGHT}, Recall={RECALL_WEIGHT}, "
                f"Prec={PRECISION_WEIGHT}, Param={PARAM_WEIGHT}\n\n")

        f.write("Metric Definitions:\n")
        f.write("  TP (True Positive): Batch step matches manual within tolerance\n")
        f.write("  FP (False Positive): Batch step doesn't match any manual step\n")
        f.write("  FN (False Negative): Manual step not detected by batch\n")
        f.write("  Precision = TP/(TP+FP): Of detected steps, how many correct?\n")
        f.write("  Recall = TP/(TP+FN): Of true steps, how many found?\n")
        f.write("  F1 = 2×P×R/(P+R): Harmonic mean of precision and recall\n\n")

        for idx, (_, row) in enumerate(results_df.iterrows()):
            f.write(f"#{idx+1} - {row['analysis_folder']}\n")
            f.write(f"   Overall Score: {row['score']}/100\n")
            f.write(f"   Classification Metrics:\n")
            f.write(f"     Precision: {row['precision']:.3f} | Recall: {row['recall']:.3f} | F1: {row['f1_score']:.3f}\n")
            f.write(f"   Contingency Table:\n")
            f.write(f"     TP={row['tp']} | FP={row['fp']} | FN={row['fn']}\n")
            if row['mean_force_dev_pN'] != 'N/A':
                f.write(f"   Parameter Accuracy:\n")
                f.write(f"     Mean Force Dev: {row['mean_force_dev_pN']} pN\n")
                f.write(f"     Mean Dist Dev: {row['mean_dist_dev_nm']} nm\n")
            f.write(f"   Parameters: {row['params_json']}\n")
            f.write(f"   Path: {row['folder_path']}\n")
            f.write("-" * 70 + "\n")

    print(f"Human-readable report saved to: {report_file}")

    # 6. Highlight top recommendations
    print("\n" + "=" * 80)
    print("TOP 3 RECOMMENDATIONS")
    print("=" * 80)

    top_k = min(3, len(results_df))
    for idx in range(top_k):
        row = results_df.iloc[idx]
        print(f"\nRank #{idx+1}: {row['analysis_folder']}")
        print(f"  Score: {row['score']}/100 | F1: {row['f1_score']}")
        print(f"  Precision: {row['precision']} | Recall: {row['recall']}")
        print(f"  TP={row['tp']} FP={row['fp']} FN={row['fn']}")
        if row['mean_force_dev_pN'] != 'N/A':
            print(f"  Avg Force Dev: {row['mean_force_dev_pN']} pN")
        print(f"  Params: {format_params(row['params_json'])}")

    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()