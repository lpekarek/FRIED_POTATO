"""FRIED POTATO -- 2024 -- Version 1"""

from scipy import signal
import numpy as np


def preprocess_RAW(Force, Distance, input_text, input_checkbox):

    if input_checkbox['aug'] == 1 and input_text['augment_factor'] != '':
        # Augment data
        new_f = []
        new_d = []
        factor = int(input_text['augment_factor'])

        for line in range(len(Force) - 1):
            f = Force[line]
            d = Distance[line]
            new_f.append(f)
            new_d.append(d)
            f_next = Force[line + 1]
            d_next = Distance[line + 1]
            delta_f = abs(f_next - f)
            delta_d = abs(d_next - d)

            for point in range(1, factor):
                mu = 0
                sigma_f = (1.5 * delta_f) / 3
                sigma_d = (1.5 * delta_d) / 3
                nf = np.random.normal(mu, sigma_f)
                nd = np.random.normal(mu, sigma_d)
                new_point_f = f + ((point + nf) / factor) * delta_f
                new_point_d = d + ((point + nd) / factor) * delta_d
                new_f.append(new_point_f)
                new_d.append(new_point_d)

        Force = np.array(new_f)
        Distance = np.array(new_d)

    if input_checkbox['prepro'] == 1:
        # Downsample
        Force_ds = Force[::input_text['downsample']]
        Distance_ds = Distance[::input_text['downsample']]

        # Filter
        b, a = signal.butter(input_text['filter_degree'], input_text['filter_cut_off'])
        filteredForce = signal.filtfilt(b, a, Force_ds)
        filteredDistance = signal.filtfilt(b, a, Distance_ds)

        Force_Distance = np.column_stack((filteredForce, filteredDistance * 1000))
        Force_Distance_um = np.column_stack((filteredForce, filteredDistance))

    else:
        Distance_nm=np.array(Distance)*1000
        Force_Distance = np.column_stack((Force, Distance_nm))
        Force_Distance_um = np.column_stack((Force, Distance))

    return Force_Distance, Force_Distance_um


# creates a dataset from min force threshold to max force value
def trim_data(FD_data, F_min):
    F_trimmed = np.array([])
    PD_trimmed = np.array([])
    F_low = np.array([])

    F_max = np.where(FD_data[:, 0] == max(FD_data[:, 0]))
    fi = F_max[0][0]

    while FD_data[fi, 0] > F_min and fi > 0:
        fi = fi - 1

    if not fi == F_max[0][0]:
        F_trimmed = FD_data[fi:F_max[0][0], 0]
        PD_trimmed = FD_data[fi:F_max[0][0], 1]
        F_low = FD_data[:fi, 0]
    elif fi == 0:
        print('Could not trim this curve, data below minimum force threshold!')

    return F_trimmed, PD_trimmed, F_low


# creates derivatives for Force and Distance of the trimmed datasets
def create_derivative(input_text, Frequency_value, F_trimmed, PD_trimmed, F_low):
    d_time = 1 / Frequency_value * input_text['downsample'] * input_text['step_d']

    x = input_text['step_d']

    derivative_list = []

    while x < len(F_trimmed):
        if PD_trimmed[0] < PD_trimmed[-1]:
            F_value = (F_trimmed[x] + F_trimmed[x - input_text['step_d']]) / 2
            PD_value = (PD_trimmed[x] + PD_trimmed[x - input_text['step_d']]) / 2
            delta_PD = PD_trimmed[x] - PD_trimmed[x - input_text['step_d']]
            delta_F = F_trimmed[x] - F_trimmed[x - input_text['step_d']]
            F_dt = delta_F / d_time
            PD_dt = delta_PD / d_time
        else:
            F_value = (F_trimmed[x] + F_trimmed[(x - input_text['step_d'])]) / 2
            PD_value = (PD_trimmed[x] + PD_trimmed[(x - input_text['step_d'])]) / 2
            delta_PD = PD_trimmed[x] - PD_trimmed[(x - input_text['step_d'])]
            delta_F = F_trimmed[x] - F_trimmed[(x - input_text['step_d'])]
            F_dt = (delta_F / d_time) * (-1)
            PD_dt = (delta_PD / d_time) * (-1)

        derivative_list.append([F_value, PD_value, F_dt, PD_dt])
        x = x + input_text['step_d']

    derivative_array = np.array(derivative_list)

    return derivative_array
