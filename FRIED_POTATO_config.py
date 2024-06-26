"""FRIED POTATO -- 2024 -- Version 1"""

"""default values for each data type"""

### Explanation - step detection parameters - depend on the input data type
# 'downsample' - Downsampling rate -- only keep each n-th data point
# 'filter_degree' - Butterworth filter degree -- for smoothing the data
# 'filter_cut_off' - Cut-off frequency of the Butterworth filter
# 'force_min' - Force threshold [pN] -- for removing data points below this threshold (non-reliable step detection at low forces)
# 'z_score_force' - Z-score threshold for force derivative -- only detect steps that exceed this threshold
# 'z_score_distance' - Z-score threshold for distance derivative -- only detect steps that exceed this threshold
# 'step d' - number of datapoints to calculate the derivative
# 'window_size' - window size (number of datapoints) for the moving median +/- Z-score filter
# 'std_difference' - for iterative data sorting
# 'frequency' - data frequency [Hz] - if available, this value will be obtained from the h5 file
# 'augment_factor' - for very low frequency data - number of data points to add to the dataset between existing data points (0 = no augmentation)

default_values_HF = {
    'downsample': '30',
    'filter_degree': '4',
    'filter_cut_off': '0.005',
    'force_min': '5',
    'z_score_force': '3',
    'z_score_distance': '3',
    'step_d': '10',
    'window_size': '800',
    'std_difference': '0.05',
    'frequency': '1000',
    'augment_factor': 0
}

default_values_LF = {
    'downsample': '1',
    'filter_degree': '2',
    'filter_cut_off': '0.5',
    'force_min': '5',
    'z_score_force': '3',
    'z_score_distance': '3',
    'step_d': '3',
    'window_size': '20',
    'std_difference': '0.05',
    'frequency': '1000',
    'augment_factor': 0
}

default_values_CSV = {
    'downsample': '1',
    'filter_degree': '1',
    'filter_cut_off': '0.01',
    'force_min': '5',
    'z_score_force': '2.5',
    'z_score_distance': '3',
    'step_d': '10',
    'window_size': '120',
    'std_difference': '0.05',
    'frequency': '20',
    'augment_factor': 10
}


### Explanation - fitting parameters - depend on the measured construct
### parameters and their maximum and minimum values allowed for fitting
# 'dsLp' - Persistence length [nm] of the double-stranded part of the construct
# 'dsLp_up' - Upper bound for the persistence length [nm] of the double-stranded part
# 'dsLp_low' - Lower bound for the persistence length [nm] of the double-stranded part
# 'ssLp' - Persistence length [nm] of the single-stranded part of the construct
# 'dsLc' - Contour length [nm] of the double-stranded part of the construct
# 'ssLc' - Contour length [nm] of the single-stranded part of the construct
# 'ssLc_up' - Upper bound for the contour length [nm] of the single-stranded part
# 'stiffness_ds' - Stiffness [pN] of the double-stranded part of the construct
# 'stiffness_ds_up' - Upper bound for the stiffness [pN] of the double-stranded part
# 'stiffness_ds_low' - Lower bound for the stiffness [pN] of the double-stranded part
# 'stiffness_ss' - Stiffness [pN] of the single-stranded part of the construct
# 'stiffness_ss_up' - Upper bound for the stiffness [pN] of the single-stranded part
# 'stiffness_ss_low' - Lower bound for the stiffness [pN] of the single-stranded part
# 'f_offset' - Force offset [pN] -- for correcting the force signal
# 'f_offset_up' - Upper bound for the force offset [pN]
# 'f_offset_low' - Lower bound for the force offset [pN]
# 'd_offset' - Distance offset [nm] -- for correcting the distance signal
# 'd_offset_up' - Upper bound for the distance offset [nm]
# 'd_offset_low' - Lower bound for the distance offset [nm]

default_values_FIT = {
    'dsLp': '40',
    'dsLp_up': '80',
    'dsLp_low': '12',
    'ssLp': '1',
    'dsLc': '1261',
    'ssLc': '1',
    'ssLc_up': '1350',
    'stiffness_ds': '500',
    'stiffness_ds_up': '600',
    'stiffness_ds_low': '400',
    'stiffness_ss': '1000',
    'stiffness_ss_up': '1500',
    'stiffness_ss_low': '300',
    'f_offset': '0',
    'f_offset_up': '0.5',
    'f_offset_low': '-0.5',
    'd_offset': '0',
    'd_offset_up': '500',
    'd_offset_low': '-500'
}

default_values_constantF = {
    'x_min_constant_f': '0',
    'x_max_constant_f': '180',
    'y_min_constant_f': '1290',
    'y_max_constant_f': '1320',
    'number_gauss': '3',
    'mean_gauss': '1295,1310,1317',
    'STD_gauss': '1,1,1',
    'amplitude_gauss': '10000,5000,10000'
}
