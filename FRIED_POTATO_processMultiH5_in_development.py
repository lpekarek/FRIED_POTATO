"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


import pandas as pd
import numpy as np

def split_H5(FD, FD_ds, input_settings, Frequency_value):
    """
    Enhanced split_H5 with Minimum Peak Force and Pause Time validation.
    
    Logic:
    1. Calculates derivative to detect direction (Pulling vs Relaxing).
    2. Uses Hysteresis to ignore transient spikes (unfolding events).
    3. Validates Cycle Transitions:
       - Must reach a Minimum Peak Force.
       - Must have a stationary "Pause" period (low derivative) for at least min_pause_time.
    
    Args:
        FD: Force-Distance array (Force, Distance)
        FD_ds: Downsampled Force-Distance array
        input_settings: Dictionary containing 'downsample_value', 'step_d'
        Frequency_value: Raw sampling frequency (Hz)
        
    Returns:
        fw_merge_keep_FD, rv_merge_keep_FD, fw_merge_keep_FDds, rv_merge_keep_FDds
    """
    
    # --- Configuration Parameters ---
    # These can be hardcoded or passed via input_settings if you prefer GUI control
    derivative_threshold_high = 3500.0  # Threshold to trigger state change (Strong signal)
    derivative_threshold_low = 1500.0   # Hysteresis threshold (Weak signal tolerance)
    
    min_peak_force = 15.0               # Minimum Force (pN) required to consider a point a "Peak"
    min_pause_time = 0.5                # Minimum time (seconds) the system must be stationary at the peak
    
    min_cycle_points = 30               # Minimum points in a segment to be valid (filters noise)
    min_cycle_displacement = 20.0       # Minimum distance change (nm) to be a valid cycle segment
    
    # --- Basic Checks ---
    if FD is None or FD_ds is None:
        raise ValueError("FD and FD_ds must be provided.")
    if len(FD) != len(FD_ds):
        raise ValueError(f"FD and FD_ds must have the same length. Got {len(FD)} vs {len(FD_ds)}.")

    # Accessors
    def row_at(arr, idx):
        return arr.iloc[idx].values if isinstance(arr, pd.DataFrame) else arr[idx]

    def ncols(arr):
        return arr.shape[1] if not isinstance(arr, pd.DataFrame) else arr.shape[1]

    # Time step calculation
    # d_time is the time step between points in the DERIVATION array
    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']
    d = Frequency_value // input_settings['downsample_value']
    d_half = int(0.5 * d)
    
    # Calculate how many points constitute the minimum pause time
    # Points = Time / Time_Per_Point
    points_per_second = 1.0 / d_time
    min_pause_points = int(min_pause_time * points_per_second)
    
    if min_pause_points < 1:
        min_pause_points = 1

    derivation_list = []
    t = 0

    # --- Build Derivative Array ---
    for i in range(d_half, len(FD) - d_half):
        PD_value = (row_at(FD, i + d_half)[1] + row_at(FD, i - d_half)[1]) / 2
        delta_PD = row_at(FD, i + d_half)[1] - row_at(FD, i - d_half)[1]
        PD_dt = delta_PD / d_time
        t = t + d_time
        derivation_list.append([t, PD_value, PD_dt])

    derivation_array = pd.DataFrame(derivation_list).to_numpy()

    # --- State Machine Variables ---
    # States: 
    # 0: Neutral/Noise
    # 1: Pulling (Forward)
    # 2: Pulling Peak (Stationary at high force)
    # 3: Relaxing (Reverse)
    # 4: Relaxing Peak (Stationary at low force)
    
    current_state = 0 
    current_direction = 0 # 1 for Pulling, -1 for Relaxing
    
    # Buffers for the current segment being built
    segment_indices = []
    segment_FD = []
    segment_FDds = []
    
    # Buffers for the "Pause" detection at the peak
    pause_indices = []
    pause_FD = []
    pause_FDds = []
    
    # Lists to store valid completed cycles
    valid_forward_cycles = []
    valid_reverse_cycles = []

    # Helper to finalize a segment
    def finalize_segment(direction, indices, fd_data, fd_ds_data, is_peak_transition=False):
        if len(indices) == 0:
            return

        # 1. Check Minimum Duration
        if len(indices) < min_cycle_points:
            return

        # 2. Check Minimum Displacement
        if len(fd_data) > 1:
            start_dist = fd_data[0][1]
            end_dist = fd_data[-1][1]
            displacement = abs(end_dist - start_dist)
            if displacement < min_cycle_displacement:
                return

        # 3. If this is a transition point (Peak), check Force and Pause criteria
        if is_peak_transition:
            # Check Max Force in the pause buffer
            if len(pause_FD) > 0:
                max_force_in_pause = np.max([row[0] for row in pause_FD])
                if max_force_in_pause < min_peak_force:
                    # Force too low, not a valid peak
                    return
                
                # Check Pause Duration
                if len(pause_indices) < min_pause_points:
                    # Not paused long enough
                    return
            else:
                # No pause detected, invalid transition
                return

        # Valid segment found
        if direction == 1:
            valid_forward_cycles.append({
                'indices': indices,
                'fd': np.array(fd_data),
                'fd_ds': np.array(fd_ds_data)
            })
        elif direction == -1:
            valid_reverse_cycles.append({
                'indices': indices,
                'fd': np.array(fd_data),
                'fd_ds': np.array(fd_ds_data)
            })

    # --- Main Loop ---
    for i in range(len(derivation_array)):
        idx = i + d_half
        deriv_val = derivation_array[i, 2]
        current_force = row_at(FD, idx)[0]
        current_dist = row_at(FD, idx)[1]
        
        # Determine Target Direction based on Hysteresis
        target_dir = 0
        if deriv_val > derivative_threshold_high:
            target_dir = 1 # Strong Pulling
        elif deriv_val < -derivative_threshold_high:
            target_dir = -1 # Strong Relaxing
        elif deriv_val > derivative_threshold_low:
            target_dir = 1 # Weak Pulling
        elif deriv_val < -derivative_threshold_low:
            target_dir = -1 # Weak Relaxing
        
        # State Transitions
        if current_state == 0:
            # Starting from noise
            if target_dir != 0:
                current_direction = target_dir
                current_state = 1 if target_dir == 1 else 3 # 1=Pulling, 3=Relaxing
                
                # Reset buffers
                segment_indices = [idx]
                segment_FD = [row_at(FD, idx)]
                segment_FDds = [row_at(FD_ds, idx)]
                pause_indices = []
                pause_FD = []
                pause_FDds = []
        
        elif current_state == 1: # Currently Pulling
            if target_dir == -1:
                # Potential Peak detected (Transition to Relaxing)
                # First, finalize the current Pulling segment IF it meets criteria
                # But wait, we need to check the PAUSE first.
                # We enter a "Checking Pause" state implicitly by collecting pause data
                
                # If we have accumulated pause data, check it now
                if len(pause_indices) >= min_pause_points:
                    finalize_segment(1, segment_indices, segment_FD, segment_FDds, is_peak_transition=True)
                
                # Start Relaxing
                current_direction = -1
                current_state = 3
                
                # Reset buffers for new segment
                segment_indices = [idx]
                segment_FD = [row_at(FD, idx)]
                segment_FDds = [row_at(FD_ds, idx)]
                pause_indices = []
                pause_FD = []
                pause_FDds = []
                
            elif target_dir == 0:
                # Noise or Stationary (Potential Pause)
                # Check if we are stationary (deriv near 0)
                if abs(deriv_val) < derivative_threshold_low:
                    # Add to pause buffer
                    pause_indices.append(idx)
                    pause_FD.append(row_at(FD, idx))
                    pause_FDds.append(row_at(FD_ds, idx))
                    
                    # Also add to current segment to maintain continuity? 
                    # Usually, we want the segment to end at the peak. 
                    # Let's add to segment but mark it as part of the peak region.
                    segment_indices.append(idx)
                    segment_FD.append(row_at(FD, idx))
                    segment_FDds.append(row_at(FD_ds, idx))
                else:
                    # Still pulling, just weak
                    segment_indices.append(idx)
                    segment_FD.append(row_at(FD, idx))
                    segment_FDds.append(row_at(FD_ds, idx))
                    # Clear pause buffer if we start moving again
                    pause_indices = []
                    pause_FD = []
                    pause_FDds = []
            else:
                # Still pulling strongly
                segment_indices.append(idx)
                segment_FD.append(row_at(FD, idx))
                segment_FDds.append(row_at(FD_ds, idx))
                pause_indices = []
                pause_FD = []
                pause_FDds = []

        elif current_state == 3: # Currently Relaxing
            if target_dir == 1:
                # Potential Peak detected (Transition to Pulling)
                if len(pause_indices) >= min_pause_points:
                    finalize_segment(-1, segment_indices, segment_FD, segment_FDds, is_peak_transition=True)
                
                current_direction = 1
                current_state = 1
                
                segment_indices = [idx]
                segment_FD = [row_at(FD, idx)]
                segment_FDds = [row_at(FD_ds, idx)]
                pause_indices = []
                pause_FD = []
                pause_FDds = []
                
            elif target_dir == 0:
                if abs(deriv_val) < derivative_threshold_low:
                    pause_indices.append(idx)
                    pause_FD.append(row_at(FD, idx))
                    pause_FDds.append(row_at(FD_ds, idx))
                    segment_indices.append(idx)
                    segment_FD.append(row_at(FD, idx))
                    segment_FDds.append(row_at(FD_ds, idx))
                else:
                    segment_indices.append(idx)
                    segment_FD.append(row_at(FD, idx))
                    segment_FDds.append(row_at(FD_ds, idx))
                    pause_indices = []
                    pause_FD = []
                    pause_FDds = []
            else:
                segment_indices.append(idx)
                segment_FD.append(row_at(FD, idx))
                segment_FDds.append(row_at(FD_ds, idx))
                pause_indices = []
                pause_FD = []
                pause_FDds = []

    # Finalize the last segment if it was a valid cycle (not just a pause)
    # We only finalize if we ended in a moving state (1 or 3) and had enough points
    if current_state in [1, 3]:
        # If we ended in a pause, we might have a valid peak but no subsequent cycle.
        # We finalize the segment leading up to the pause.
        if len(segment_indices) >= min_cycle_points:
             # Check displacement
            if len(segment_FD) > 1:
                disp = abs(segment_FD[-1][1] - segment_FD[0][1])
                if disp >= min_cycle_displacement:
                     if current_direction == 1:
                        valid_forward_cycles.append({'indices': segment_indices, 'fd': np.array(segment_FD), 'fd_ds': np.array(segment_FDds)})
                     else:
                        valid_reverse_cycles.append({'indices': segment_indices, 'fd': np.array(segment_FD), 'fd_ds': np.array(segment_FDds)})

    # --- Reconstruct Output Arrays ---
    fd_cols = ncols(FD)
    ds_cols = ncols(FD_ds)

    fw_list = [seg['fd'] for seg in valid_forward_cycles]
    rv_list = [seg['fd'] for seg in valid_reverse_cycles]
    fw_ds_list = [seg['fd_ds'] for seg in valid_forward_cycles]
    rv_ds_list = [seg['fd_ds'] for seg in valid_reverse_cycles]

    # Handle empty cases
    if len(fw_list) == 0: fw_list = [np.empty((0, fd_cols))]
    if len(rv_list) == 0: rv_list = [np.empty((0, fd_cols))]
    if len(fw_ds_list) == 0: fw_ds_list = [np.empty((0, ds_cols))]
    if len(rv_ds_list) == 0: rv_ds_list = [np.empty((0, ds_cols))]

    fw_merge_keep_FD = np.array(fw_list, dtype=object)
    rv_merge_keep_FD = np.array(rv_list, dtype=object)
    fw_merge_keep_FDds = np.array(fw_ds_list, dtype=object)
    rv_merge_keep_FDds = np.array(rv_ds_list, dtype=object)

    print(f"Split Results: {len(fw_list)} Forward, {len(rv_list)} Reverse cycles.")
    print(f"Criteria: Min Peak Force={min_peak_force}pN, Min Pause={min_pause_time}s")

    return fw_merge_keep_FD, rv_merge_keep_FD, fw_merge_keep_FDds, rv_merge_keep_FDds


#version from 04-05-2026
def split_H5_not_so_old(FD, FD_ds, input_settings, Frequency_value):
    """
    Split FD into forward/reverse segments based on time derivative,
    and split FD_ds into the exact same segments using the same indices.

    Returns:
        fw_merge_keep_FD, rv_merge_keep_FD, fw_merge_keep_FDds, rv_merge_keep_FDds
        (each is an array(dtype=object) of segment arrays)
    """
    derivative_threshold = 2500
    # --- Basic checks ---
    if FD is None or FD_ds is None:
        raise ValueError("FD and FD_ds must be provided.")
    if len(FD) != len(FD_ds):
        raise ValueError(f"FD and FD_ds must have the same length. Got {len(FD)} vs {len(FD_ds)}.")

    # Accessors to handle DataFrame vs ndarray seamlessly
    def row_at(arr, idx):
        return arr.iloc[idx].values if isinstance(arr, pd.DataFrame) else arr[idx]

    def ncols(arr):
        return arr.shape[1] if not isinstance(arr, pd.DataFrame) else arr.shape[1]

    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']
    d = Frequency_value // input_settings['downsample_value']
    d_half = int(0.5 * d)
    derivation_list = []
    t = 0

    # --- Build derivative from FD only ---
    for i in range(d_half, len(FD) - d_half):
        # Column 1 is assumed to be PD (same as your code)
        PD_value = (row_at(FD, i + d_half)[1] + row_at(FD, i - d_half)[1]) / 2
        delta_PD = row_at(FD, i + d_half)[1] - row_at(FD, i - d_half)[1]
        PD_dt = delta_PD / d_time
        t = t + d_time
        derivation_list.append([t, PD_value, PD_dt])

    derivation_array = pd.DataFrame(derivation_list).to_numpy()

    forward_FD = []
    reverse_FD = []
    forward_FDds = []
    reverse_FDds = []

    n = []
    x_fw = []
    x_rv = []
    x_total = []

    # --- Classify points & collect the same rows from FD_ds by index ---
    for i in range(len(derivation_array)):
        idx = i + d_half
        if -derivative_threshold < derivation_array[i, 2] < derivative_threshold:
            n.append(idx)
        elif derivation_array[i, 2] > derivative_threshold:
            x_fw.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            forward_FD.append(row_at(FD, idx))
            forward_FDds.append(row_at(FD_ds, idx))
        elif derivation_array[i, 2] < -derivative_threshold:
            x_rv.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            reverse_FD.append(row_at(FD, idx))
            reverse_FDds.append(row_at(FD_ds, idx))

    unique_fw = np.unique(x_fw)
    unique_rv = np.unique(x_rv)
    unique_total = np.unique(x_total)
    print('Fw:', len(unique_fw), ', Rev:', len(unique_rv), ', Together:', len(unique_total))

    # --- Stack (handle empty) ---
    fd_cols = ncols(FD)
    ds_cols = ncols(FD_ds)

    forward_FD = np.vstack(forward_FD) if len(forward_FD) > 0 else np.empty((0, fd_cols))
    reverse_FD = np.vstack(reverse_FD) if len(reverse_FD) > 0 else np.empty((0, fd_cols))

    forward_FDds = np.vstack(forward_FDds) if len(forward_FDds) > 0 else np.empty((0, ds_cols))
    reverse_FDds = np.vstack(reverse_FDds) if len(reverse_FDds) > 0 else np.empty((0, ds_cols))

    # --- Attach the index labels (same x_fw/x_rv) for both FD and FD_ds ---
    fw_merge_FD   = np.column_stack((forward_FD,  x_fw)) if len(forward_FD)   > 0 else np.empty((0, fd_cols + 1))
    rv_merge_FD   = np.column_stack((reverse_FD,  x_rv)) if len(reverse_FD)   > 0 else np.empty((0, fd_cols + 1))
    fw_merge_FDds = np.column_stack((forward_FDds, x_fw)) if len(forward_FDds) > 0 else np.empty((0, ds_cols + 1))
    rv_merge_FDds = np.column_stack((reverse_FDds, x_rv)) if len(reverse_FDds) > 0 else np.empty((0, ds_cols + 1))

    # --- Split into contiguous blocks by changes in the last column ---
    # (uses the exact same boundaries for FD and FD_ds)
    def split_by_label(arr):
        if arr.shape[0] == 0:
            return []
        cuts = np.where(np.diff(arr[:, -1]))[0] + 1
        return np.split(arr, cuts)

    fw_merge_FD   = split_by_label(fw_merge_FD)
    rv_merge_FD   = split_by_label(rv_merge_FD)
    fw_merge_FDds = split_by_label(fw_merge_FDds)
    rv_merge_FDds = split_by_label(rv_merge_FDds)

    # --- Keep only segments longer than half of average segment length (based on FD) ---
    arr_length = sum(len(i) for i in fw_merge_FD + rv_merge_FD)
    count = len(fw_merge_FD) + len(rv_merge_FD)
    arr_average = arr_length / count if count > 0 else 0

    def keep_long_enough(blocks):
        return [b for b in blocks if len(b) >= 0.5 * arr_average]

    fw_keep_FD   = keep_long_enough(fw_merge_FD)
    rv_keep_FD   = keep_long_enough(rv_merge_FD)
    # Apply the same length threshold to FD_ds
    fw_keep_FDds = keep_long_enough(fw_merge_FDds)
    rv_keep_FDds = keep_long_enough(rv_merge_FDds)

    # --- Return without the label column for the user data, keeping arrays as ragged (dtype=object) ---
    def strip_label(blocks):
        return [b[:, :-1] if b.shape[1] > 0 else b for b in blocks]

    return np.array(strip_label(fw_keep_FD),   dtype=object), np.array(strip_label(rv_keep_FD),   dtype=object), np.array(strip_label(fw_keep_FDds), dtype=object),np.array(strip_label(rv_keep_FDds), dtype=object),




#original version of split_H5 12-08-2025
def split_H5_old(FD, FD_ds, input_settings, Frequency_value):
    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']
    d = Frequency_value // input_settings['downsample_value']
    d_half = int(0.5 * d)
    derivation_list = []
    t = 0

    for i in range(d_half, len(FD) - d_half):
        PD_value = (FD[i + d_half, 1] + FD[i - d_half, 1]) / 2
        delta_PD = FD[i + d_half, 1] - FD[i - d_half, 1]
        PD_dt = delta_PD / d_time
        t = t + d_time
        derivation_list.append([t, PD_value, PD_dt])

    derivation_array = pd.DataFrame(derivation_list).to_numpy()

    forward = []
    reverse = []
    n = []
    x_fw = []
    x_rv = []
    x_total = []

    for i in range(len(derivation_array)):
        if -2500 < derivation_array[i, 2] < 2500:
            n.append(i + d_half)
        elif derivation_array[i, 2] > 2500:
            x_fw.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            forward.append(FD[i + d_half].values if isinstance(FD, pd.DataFrame) else FD[i + d_half])
        elif derivation_array[i, 2] < -2500:
            x_rv.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            reverse.append(FD[i + d_half].values if isinstance(FD, pd.DataFrame) else FD[i + d_half])

    unique_fw = np.unique(x_fw)
    unique_rv = np.unique(x_rv)
    unique_total = np.unique(x_total)

    print('Fw:', len(unique_fw), ', Rev:', len(unique_rv), ', Together:', len(unique_total))

    # Ensure arrays are not empty before stacking
    forward = np.vstack(forward) if len(forward) > 0 else np.empty((0, FD.shape[1]))
    reverse = np.vstack(reverse) if len(reverse) > 0 else np.empty((0, FD.shape[1]))

    print(forward)

    fw_merge = np.column_stack((forward, x_fw)) if len(forward) > 0 else np.empty((0, FD.shape[1] + 1))
    rv_merge = np.column_stack((reverse, x_rv)) if len(reverse) > 0 else np.empty((0, FD.shape[1] + 1))

    fw_merge = np.split(fw_merge, np.where(np.diff(fw_merge[:, -1]))[0] + 1)
    rv_merge = np.split(rv_merge, np.where(np.diff(rv_merge[:, -1]))[0] + 1)

    # Remove arrays that are too short
    arr_length = sum(len(i) for i in fw_merge + rv_merge)
    count = len(fw_merge) + len(rv_merge)
    arr_average = arr_length / count if count > 0 else 0

    fw_merge_keep = [i for i in fw_merge if len(i) >= 0.5 * arr_average]
    rv_merge_keep = [i for i in rv_merge if len(i) >= 0.5 * arr_average]

    return np.array(fw_merge_keep, dtype=object), np.array(rv_merge_keep, dtype=object)


"""def split_H5(FD, input_settings, Frequency_value):
    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']
    d = Frequency_value // input_settings['downsample_value']
    d_half = int(0.5 * d)
    derivation_list = []
    t = 0

    for i in range(d_half, len(FD) - d_half):
        PD_value = (FD[i + d_half, 1] + FD[i - d_half, 1]) / 2
        delta_PD = FD[i + d_half, 1] - FD[i - d_half, 1]
        PD_dt = delta_PD / d_time
        t = t + d_time
        derivation_list.append([t, PD_value, PD_dt])

    derivation_array = pd.DataFrame(derivation_list)
    derivation_array = derivation_array.to_numpy()

    forward = []
    reverse = []
    n = []
    x = d_half
    x_fw = []
    x_rv = []
    x_total = []

    for i in range(len(derivation_array)):
        if derivation_array[i, 2] > -2500 and derivation_array[i, 2] < 2500:
            n.append(i + d_half)
            x = len(n) + d_half
        elif derivation_array[i, 2] > 2500:
            x_fw.append(x)
            x_total.append(x)
            forward.append(FD[i + d_half])
        elif derivation_array[i, 2] < -2500:
            x_rv.append(x)a
            x_total.append(x)
            reverse.append(FD[i + d_half])

    unique_fw = np.unique(x_fw)
    unique_rv = np.unique(x_rv)
    unique_total = np.unique(x_total)

    print('Fw:', len(unique_fw), ', Rev:', len(unique_rv), ', Together:', len(unique_total))

    forward = np.vstack(forward)
    reverse = np.vstack(reverse)
    print(forward)

    fw_merge = np.column_stack((forward, x_fw))
    rv_merge = np.column_stack((reverse, x_rv))

    fw_merge = np.split(fw_merge, np.where(np.diff(fw_merge[:, 2]))[0] + 1)
    rv_merge = np.split(rv_merge, np.where(np.diff(rv_merge[:, 2]))[0] + 1)
    
    ############################# remove arrays that are way below average length
    arr_length = 0
    count = 0
    for i in fw_merge:
        arr_length += len(i)
        count += 1
    for i in rv_merge:
        arr_length += len(i)
        count += 1
    arr_average = arr_length / count
    try: 
        for i in fw_merge:
            if len(i) < 0.5 * arr_average:
                fw_merge_keep = np.delete(fw_merge, i)
        for i in rv_merge:
            if len(i) < 0.5 * arr_average:
                rv_merge_keep = np.delete(rv_merge, i)
    except Exception as e: 
        print(e)

    return fw_merge_keep, rv_merge_keep"""
