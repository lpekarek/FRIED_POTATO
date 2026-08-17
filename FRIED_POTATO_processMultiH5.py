"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""

import pandas as pd
import numpy as np
import h5py

def check_for_trap_position(file_path):
    """
    Helper function to check if 'Trap Position 1x' exists in the H5 file.
    Returns True if present, False otherwise.
    """
    try:
        with h5py.File(file_path, "r") as f:
            # Common paths for Lumicks C-Trap
            paths_to_check = [
                "Trap position/1X",
                "Trap position/2X",
                "Distance/Piezo Distance" # Alternative naming
            ]
            
            for path in paths_to_check:
                if path in f:
                    return True
        return False
    except Exception:
        return False

def split_H5(FD, FD_ds, input_settings, Frequency_value, use_trap_pos=False, trap_data=None, trap_ds_data=None):
    """
    Split FD into forward/reverse segments based on time derivative.
    
    ENHANCED LOGIC:
    - If use_trap_pos=True: Uses Trap Position data for derivative calculation (avoids unfolding noise)
    - If use_trap_pos=False: Falls back to Distance data (legacy behavior for older files)
    
    Args:
        FD: Force-Distance array (Force, Distance) - Used for output data
        FD_ds: Downsampled Force-Distance array
        input_settings: Dictionary containing 'downsample_value', 'step_d'
        Frequency_value: Raw sampling frequency (Hz)
        use_trap_pos: Boolean flag indicating if Trap Position data is available
        trap_data: Array of Trap Position values (same length as FD, column 1 = position)
        trap_ds_data: Downsampled Trap Position values
        
    Returns:
        fw_merge_keep_FD, rv_merge_keep_FD, fw_merge_keep_FDds, rv_merge_keep_FDds
    """
    derivative_threshold = 2500
    
    # --- Basic Checks ---
    if FD is None or FD_ds is None:
        raise ValueError("FD and FD_ds must be provided.")
    if len(FD) != len(FD_ds):
        raise ValueError(f"FD and FD_ds must have the same length. Got {len(FD)} vs {len(FD_ds)}.")
    
    # Check if trap_data is valid when requested
    if use_trap_pos:
        if trap_data is None or len(trap_data) != len(FD):
            print("WARNING: use_trap_pos=True but trap_data missing/mismatched. Falling back to Distance.")
            use_trap_pos = False

    # Accessors
    def row_at(arr, idx):
        return arr.iloc[idx].values if isinstance(arr, pd.DataFrame) else arr[idx]

    def ncols(arr):
        return arr.shape[1] if not isinstance(arr, pd.DataFrame) else arr.shape[1]

    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']
    d = Frequency_value // input_settings['downsample_value']
    d_half = int(0.5 * d)
    
    derivation_list = []
    t = 0

    # --- Build Derivative from TRAP position if available, else DISTANCE ---
    source_data = trap_data if use_trap_pos else FD
    
    for i in range(d_half, len(source_data) - d_half):
        coord_value = (row_at(source_data, i + d_half)[1] + row_at(source_data, i - d_half)[1]) / 2
        delta_coord = row_at(source_data, i + d_half)[1] - row_at(source_data, i - d_half)[1]
        coord_dt = delta_coord / d_time
        t = t + d_time
        derivation_list.append([t, coord_value, coord_dt])

    derivation_array = pd.DataFrame(derivation_list).to_numpy()

    forward_FD = []
    reverse_FD = []
    forward_FDds = []
    reverse_FDds = []

    n = []
    x_fw = []
    x_rv = []
    x_total = []

    # --- Classify points & collect rows from FD_ds by index ---
    for i in range(len(derivation_array)):
        idx = i + d_half
        deriv_val = derivation_array[i, 2]
        
        if -derivative_threshold < deriv_val < derivative_threshold:
            n.append(idx)
        elif deriv_val > derivative_threshold:
            x_fw.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            forward_FD.append(row_at(FD, idx))
            forward_FDds.append(row_at(FD_ds, idx))
        elif deriv_val < -derivative_threshold:
            x_rv.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            reverse_FD.append(row_at(FD, idx))
            reverse_FDds.append(row_at(FD_ds, idx))

    unique_fw = np.unique(x_fw)
    unique_rv = np.unique(x_rv)
    unique_total = np.unique(x_total)
    
    mode_str = "[Trap Mode]" if use_trap_pos else "[Distance Mode]"
    print(f'{mode_str} Fw: {len(unique_fw)}, Rev: {len(unique_rv)}, Together: {len(unique_total)}')

    # --- Stack (handle empty) ---
    fd_cols = ncols(FD)
    ds_cols = ncols(FD_ds)

    forward_FD = np.vstack(forward_FD) if len(forward_FD) > 0 else np.empty((0, fd_cols))
    reverse_FD = np.vstack(reverse_FD) if len(reverse_FD) > 0 else np.empty((0, fd_cols))

    forward_FDds = np.vstack(forward_FDds) if len(forward_FDds) > 0 else np.empty((0, ds_cols))
    reverse_FDds = np.vstack(reverse_FDds) if len(reverse_FDds) > 0 else np.empty((0, ds_cols))

    # --- Attach index labels ---
    fw_merge_FD   = np.column_stack((forward_FD,  x_fw)) if len(forward_FD)   > 0 else np.empty((0, fd_cols + 1))
    rv_merge_FD   = np.column_stack((reverse_FD,  x_rv)) if len(reverse_FD)   > 0 else np.empty((0, fd_cols + 1))
    fw_merge_FDds = np.column_stack((forward_FDds, x_fw)) if len(forward_FDds) > 0 else np.empty((0, ds_cols + 1))
    rv_merge_FDds = np.column_stack((reverse_FDds, x_rv)) if len(reverse_FDds) > 0 else np.empty((0, ds_cols + 1))

    # --- Split into contiguous blocks by changes in last column ---
    def split_by_label(arr):
        if arr.shape[0] == 0:
            return []
        cuts = np.where(np.diff(arr[:, -1]))[0] + 1
        return np.split(arr, cuts)

    fw_merge_FD   = split_by_label(fw_merge_FD)
    rv_merge_FD   = split_by_label(rv_merge_FD)
    fw_merge_FDds = split_by_label(fw_merge_FDds)
    rv_merge_FDds = split_by_label(rv_merge_FDds)

    # --- Keep only segments longer than half of average segment length ---
    arr_length = sum(len(i) for i in fw_merge_FD + rv_merge_FD)
    count = len(fw_merge_FD) + len(rv_merge_FD)
    arr_average = arr_length / count if count > 0 else 0

    def keep_long_enough(blocks):
        return [b for b in blocks if len(b) >= 0.5 * arr_average]

    fw_keep_FD   = keep_long_enough(fw_merge_FD)
    rv_keep_FD   = keep_long_enough(rv_merge_FD)
    fw_keep_FDds = keep_long_enough(fw_merge_FDds)
    rv_keep_FDds = keep_long_enough(rv_merge_FDds)

    # --- Return without the label column (ragged arrays) ---
    def strip_label(blocks):
        return [b[:, :-1] if b.shape[1] > 0 else b for b in blocks]

    return (np.array(strip_label(fw_keep_FD), dtype=object), 
            np.array(strip_label(rv_keep_FD), dtype=object), 
            np.array(strip_label(fw_keep_FDds), dtype=object),
            np.array(strip_label(rv_keep_FDds), dtype=object))