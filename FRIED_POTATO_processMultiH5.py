"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


import pandas as pd
import numpy as np


def split_H5(FD, input_settings, Frequency_value):
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
