"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


import pandas as pd
import numpy as np

def split_H5(FD, FD_ds, input_settings, Frequency_value):
    """
    Split FD into forward/reverse segments based on time derivative,
    and split FD_ds into the exact same segments using the same indices.

    Returns:
        fw_merge_keep_FD, rv_merge_keep_FD, fw_merge_keep_FDds, rv_merge_keep_FDds
        (each is an array(dtype=object) of segment arrays)
    """
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
        if -2500 < derivation_array[i, 2] < 2500:
            n.append(idx)
        elif derivation_array[i, 2] > 2500:
            x_fw.append(len(n) + d_half)
            x_total.append(len(n) + d_half)
            forward_FD.append(row_at(FD, idx))
            forward_FDds.append(row_at(FD_ds, idx))
        elif derivation_array[i, 2] < -2500:
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
            x_rv.append(x)
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
