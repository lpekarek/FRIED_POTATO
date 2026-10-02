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

def calculate_metrics(manual_df, batch_df, 
                      force_tol=0.3, 
                      dist_tol=1.0):
    """
    Calculate detailed match metrics between manual and batch.
    Returns dict with all components for transparency.
    """
    metrics_per_curve = {}
    curves = manual_df['filename'].unique()
    
    for curve in curves:
        m_rows = manual_df[manual_df['filename'] == curve].sort_values('step number')
        b_rows = batch_df[batch_df['filename'] == curve].sort_values('step number')
        
        n_manual = len(m_rows)
        n_batch = len(b_rows)
        
        # --- STEP COUNT COMPARISON ---
        n_missed = max(0, n_manual - n_batch)  # False negatives
        n_spurious = max(0, n_batch - n_manual)  # False positives
        
        # --- PARAMETER MATCHING FOR OVERLAPPING STEPS ---
        force_devs = []
        dist_devs = []
        matched_steps = 0
        
        min_len = min(n_manual, n_batch)
        for i in range(min_len):
            m_row = m_rows.iloc[i]
            b_row = b_rows.iloc[i]
            
            m_force = m_row.get('mean force [pn]', np.nan)
            b_force = b_row.get('fc', np.nan)
            
            m_dist = m_row.get('step length [nm]', np.nan)
            b_dist = b_row.get('step length', np.nan)
            
            if pd.notna(m_force) and pd.notna(b_force):
                force_devs.append(abs(m_force - b_force))
            
            if pd.notna(m_dist) and pd.notna(b_dist):
                dist_devs.append(abs(m_dist - b_dist))
            
            if (pd.notna(m_force) and pd.notna(b_force) and 
                pd.notna(m_dist) and pd.notna(b_dist)):
                matched_steps += 1
        
        metrics_per_curve[curve] = {
            'n_manual': n_manual,
            'n_batch': n_batch,
            'n_missed': n_missed,
            'n_spurious': n_spurious,
            'mean_force_dev_pN': np.mean(force_devs) if force_devs else np.nan,
            'mean_dist_dev_nm': np.mean(dist_devs) if dist_devs else np.nan,
            'matched_steps': matched_steps
        }
    
    return metrics_per_curve

def compute_score(metrics_per_curve,
                  fn_weight=1.0,      # Weight for false negatives (missed steps)
                  fp_weight=0.5,       # Weight for false positives (spurious steps)
                  param_weight=0.1,    # Weight for parameter deviation
                  force_tol=0.3,
                  dist_tol=1.0):
    """
    Compute overall score based on biological priorities.
    Higher is better. Score is normalized to 0-100 range.
    
    Parameters:
    -----------
    fn_weight : float
        Penalty multiplier for missed steps (false negatives)
    fp_weight : float  
        Penalty multiplier for spurious steps (false positives)
    param_weight : float
        Penalty multiplier for force/distance deviations beyond tolerance
    """
    total_penalty = 0.0
    n_curves = len(metrics_per_curve)
    
    # Collect deviations for normalization
    all_force_devs = []
    all_dist_devs = []
    
    for curve, m in metrics_per_curve.items():
        # Step count penalties (primary concern)
        total_penalty += m['n_missed'] * fn_weight * 10  # 10 pts per missed step
        total_penalty += m['n_spurious'] * fp_weight * 5  # 5 pts per spurious step (half penalty)
        
        # Parameter deviation penalties (only if within reasonable range)
        if pd.notna(m['mean_force_dev_pN']):
            excess_force = max(0, m['mean_force_dev_pN'] - force_tol)
            all_force_devs.append(excess_force)
        
        if pd.notna(m['mean_dist_dev_nm']):
            excess_dist = max(0, m['mean_dist_dev_nm'] - dist_tol)
            all_dist_devs.append(excess_dist)
    
    # Add parameter deviation penalty (normalized)
    if all_force_devs:
        avg_force_excess = np.mean(all_force_devs)
        total_penalty += avg_force_excess * param_weight * 100
    
    if all_dist_devs:
        avg_dist_excess = np.mean(all_dist_devs)
        total_penalty += avg_dist_excess * param_weight * 100
    
    # Convert to 0-100 scale (100 = perfect match)
    max_possible_penalty = n_curves * (max([m['n_manual'] for m in metrics_per_curve.values()]) * fn_weight * 10)
    
    if max_possible_penalty == 0:
        return 100.0
    
    score = 100.0 * (1 - total_penalty / max_possible_penalty)
    return max(0.0, min(100.0, score))

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
    
    # Scoring tolerances (tune based on your data quality)
    FORCE_TOL_PN = 0.3      # +/- pN tolerance for "close enough"
    DIST_TOL_NM = 1.0       # +/- nm tolerance for "close enough"
    
    # Biological prioritization (tune these!)
    FN_WEIGHT = 1.0         # Missed step penalty multiplier (higher = stricter on missing steps)
    FP_WEIGHT = 0.5         # Spurious step penalty multiplier (lower = more forgiving on extra steps)
    PARAM_WEIGHT = 0.1      # Parameter deviation importance (0 = ignore, 1 = equal to step count)
    
    print("=" * 80)
    print("FRIED POTATO OPTIMIZER -- Manual Analysis Match Focus")
    print("=" * 80)
    print(f"\nConfiguration:")
    print(f"  Root path: {ROOT_PATH}")
    print(f"  Force tolerance: ±{FORCE_TOL_PN} pN")
    print(f"  Distance tolerance: ±{DIST_TOL_NM} nm")
    print(f"  FN weight (missed steps): {FN_WEIGHT}")
    print(f"  FP weight (spurious steps): {FP_WEIGHT}")
    print(f"  Param weight (precision): {PARAM_WEIGHT}")
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
    
    # 3. Score each run with detailed diagnostics
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
        
        # Calculate metrics
        metrics = calculate_metrics(manual_df, batch_df, 
                                    force_tol=FORCE_TOL_PN, 
                                    dist_tol=DIST_TOL_NM)
        
        # Compute overall score
        score = compute_score(metrics,
                              fn_weight=FN_WEIGHT,
                              fp_weight=FP_WEIGHT,
                              param_weight=PARAM_WEIGHT,
                              force_tol=FORCE_TOL_PN,
                              dist_tol=DIST_TOL_NM)
        
        # Aggregate curve-level metrics
        total_missed = sum(m['n_missed'] for m in metrics.values())
        total_spurious = sum(m['n_spurious'] for m in metrics.values())
        mean_force_dev = np.nanmean([m['mean_force_dev_pN'] for m in metrics.values() 
                                     if pd.notna(m['mean_force_dev_pN'])])
        mean_dist_dev = np.nanmean([m['mean_dist_dev_nm'] for m in metrics.values()
                                    if pd.notna(m['mean_dist_dev_nm'])])
        
        # Count curves with perfect step count match
        n_exact_match = sum(1 for m in metrics.values() if m['n_manual'] == m['n_batch'])
        
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
            'curves_total': n_manual_curves,
            'curves_exact_match': n_exact_match,
            'steps_missed_total': total_missed,
            'steps_spurious_total': total_spurious,
            'mean_force_dev_pN': round(mean_force_dev, 2) if not np.isnan(mean_force_dev) else 'N/A',
            'mean_dist_dev_nm': round(mean_dist_dev, 2) if not np.isnan(mean_dist_dev) else 'N/A',
            'params_json': str(display_params)
        })
        
        print(f"  → Score: {score:.1f}/100 | Exact matches: {n_exact_match}/{n_manual_curves} | "
              f"Missed: {total_missed} | Spurious: {total_spurious}")
    
    # 4. Rank and display results
    if not results:
        print("\nNo valid results to rank.")
        return
    
    results_df = pd.DataFrame(results).sort_values(by='score', ascending=False)
    
    print("\n" + "=" * 80)
    print("RANKING BY MANUAL ANALYSIS MATCH QUALITY")
    print("=" * 80)
    
    # Create readable parameter summary
    def format_params(param_str):
        try:
            d = json.loads(param_str.replace("'", '"'))
            return ", ".join([f"{k}={v}" for k, v in list(d.items())[:4]])
        except:
            return param_str[:50] + "..." if len(param_str) > 50 else param_str
    
    summary_df = results_df.copy()
    summary_df['params_summary'] = summary_df['params_json'].apply(format_params)
    
    cols_to_show = ['analysis_folder', 'score', 'curves_exact_match', 
                    'steps_missed_total', 'steps_spurious_total', 'params_summary']
    print(summary_df[cols_to_show].to_string(index=False))
    
    # 5. Save results
    output_file = Path(ROOT_PATH) / "parameter_optimization_ranking.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nDetailed rankings saved to: {output_file}")
    
    # Human-readable report
    report_file = Path(ROOT_PATH) / "parameter_optimization_report.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("FRIED POTATO Parameter Optimization Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Manual Analysis: {manual_file}\n")
        f.write(f"Total Curves: {n_manual_curves}\n")
        f.write(f"Total Manual Steps: {n_manual_steps}\n")
        f.write(f"Force Tolerance: ±{FORCE_TOL_PN} pN\n")
        f.write(f"Distance Tolerance: ±{DIST_TOL_NM} nm\n\n")
        f.write(f"Scoring Weights: FN={FN_WEIGHT}, FP={FP_WEIGHT}, Param={PARAM_WEIGHT}\n\n")
        
        for idx, row in results_df.iterrows():
            f.write(f"#{idx+1} - {row['analysis_folder']}\n")
            f.write(f"   Overall Score: {row['score']}/100\n")
            f.write(f"   Curves with Exact Step Count Match: {row['curves_exact_match']}/{row['curves_total']}\n")
            f.write(f"   Total Steps Missed (FN): {row['steps_missed_total']}\n")
            f.write(f"   Total Spurious Steps (FP): {row['steps_spurious_total']}\n")
            if row['mean_force_dev_pN'] != 'N/A':
                f.write(f"   Mean Force Deviation: {row['mean_force_dev_pN']} pN\n")
                f.write(f"   Mean Distance Deviation: {row['mean_dist_dev_nm']} nm\n")
            f.write(f"   Parameters: {row['params_json']}\n")
            f.write(f"   Path: {row['folder_path']}\n")
            f.write("-" * 60 + "\n")
    
    print(f"Human-readable report saved to: {report_file}")
    
    # 6. Highlight top recommendations
    print("\n" + "=" * 80)
    print("TOP RECOMMENDATIONS")
    print("=" * 80)
    
    top_k = min(3, len(results_df))
    for idx in range(top_k):
        row = results_df.iloc[idx]
        print(f"\nRank #{idx+1}: {row['analysis_folder']}")
        print(f"  Score: {row['score']}/100")
        print(f"  Perfect Step Count Matches: {row['curves_exact_match']}/{row['curves_total']}")
        if row['mean_force_dev_pN'] != 'N/A':
            print(f"  Avg Force Dev: {row['mean_force_dev_pN']} pN | Avg Dist Dev: {row['mean_dist_dev_nm']} nm")
        print(f"  Params: {format_params(row['params_json'])}")
    
    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()