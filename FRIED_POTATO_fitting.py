"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


import matplotlib.pyplot as plt
import lumicks.pylake as lk
import numpy as np
import pandas as pd
from scipy.integrate import simps
from matplotlib.lines import Line2D
import matplotlib.patches as patches

"""define the functions used for fitting"""


# if the step start/end are selected manually by TOMATO,
# the nearest point on the curve is selected in the array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return idx


def fitting_ds(filename_i, input_settings, export_data, input_fitting, i_start, Force_Distance, derivative_array, F_low, TOMATO_param):
    global model_ds, fit_ds
    global ds_fit_dict
    global f_fitting_region_ds, d_fitting_region_ds
    global export_fit_ds
    global fitting_model


    if TOMATO_param == 0:
        start_step1 = np.where(derivative_array[:, 1] == i_start)
        start_step1 = start_step1[0][0]
        f_fitting_region_ds = Force_Distance[0:start_step1 * input_settings['step_d'] + len(F_low), 0]
        d_fitting_region_ds = Force_Distance[0:start_step1 * input_settings['step_d'] + len(F_low), 1]
    elif TOMATO_param == 1:
        start_step1 = find_nearest(Force_Distance[:, 1], i_start)
        f_fitting_region_ds = Force_Distance[0:start_step1, 0]
        d_fitting_region_ds = Force_Distance[0:start_step1, 1]

    #model_ds = lk.inverted_odijk("ds_part").subtract_independent_offset() + lk.force_offset("ds_part") #original version in POTATO 04-06-2024

    model_ds = lk.ewlc_odijk_force("ds_part").subtract_independent_offset() + lk.force_offset("ds_part") # updated version according to Lumicks changelog v0.13.2 | 2022-11-15¶ 
    fit_ds = lk.FdFit(model_ds)

    fit_ds.add_data("Double stranded", f_fitting_region_ds, d_fitting_region_ds)
    
    # Persistance length bounds
    fit_ds["ds_part/Lp"].value = input_fitting['lp_ds']
    fit_ds["ds_part/Lp"].upper_bound = input_fitting['lp_ds_up']
    fit_ds["ds_part/Lp"].lower_bound = input_fitting['lp_ds_low']

    # Force shift bounds
    fit_ds["ds_part/f_offset"].value = input_fitting['offset_f']
    fit_ds["ds_part/f_offset"].upper_bound = input_fitting['offset_f_up']
    fit_ds["ds_part/f_offset"].lower_bound = input_fitting['offset_f_low']

    # distance shift bounds
    fit_ds["ds_part/d_offset"].value = input_fitting['offset_d']
    fit_ds["ds_part/d_offset"].upper_bound = input_fitting['offset_d_up']
    fit_ds["ds_part/d_offset"].lower_bound = input_fitting['offset_d_low']
    fit_ds["ds_part/d_offset"].unit='nm'
    # stiffnes
    fit_ds["ds_part/St"].value = input_fitting['ds_stiff']
    fit_ds["ds_part/St"].upper_bound = input_fitting['ds_stiff_up']
    fit_ds["ds_part/St"].lower_bound = input_fitting['ds_stiff_low']

    # contour length
    Lc_initial_guess = input_fitting['lc_ds']  # nm
    Lc_range = 5
    fit_ds["ds_part/Lc"].upper_bound = Lc_initial_guess + Lc_range
    fit_ds["ds_part/Lc"].lower_bound = Lc_initial_guess - Lc_range
    fit_ds["ds_part/Lc"].value = Lc_initial_guess
    fit_ds["ds_part/Lc"].unit = 'nm'

    fit_ds.fit()
    fit_qual = fit_ds.log_likelihood()
    print(fit_ds.params)

    # calculate the integral until the first unfolding step
    # used to calculate the work done by the machine
    distance_integral = np.arange(min(Force_Distance[:, 1]), i_start)
    ds_integral = model_ds(distance_integral, fit_ds.params)
    area_ds = simps(ds_integral)
    print("area_ds = " + str(area_ds))

    # export the fitting parameters
    ds_fit_dict = {
        'filename': filename_i,
        'model': 'WLC',
        'model_ds': model_ds,
        'model_ss': fit_ds,
        'fit_model': fit_ds,
        'log_likelihood': fit_qual,
        'Lc_ds': fit_ds["ds_part/Lc"].value,
        'Lp_ds': fit_ds["ds_part/Lp"].value,
        'Lp_ds_stderr': fit_ds["ds_part/Lp"].stderr,
        'St_ds': fit_ds["ds_part/St"].value,
        'f_offset': fit_ds["ds_part/f_offset"].value,
        'd_offset': fit_ds["ds_part/d_offset"].value,
        'Lc_ss':0
    }
    print(ds_fit_dict)
    return ds_fit_dict, area_ds, start_step1

def fitting_FU(filename_i, input_settings, export_data, input_fitting, first_step_start,i_start,  i_end, Force_Distance, derivative_array, F_low, TOMATO_param):
    global model_ss
    global ss_fit_dict

    global model_ds, fit_ds
    global ds_fit_dict
    global f_fitting_region_ds, d_fitting_region_ds
    global export_fit_ds
    global fitting_model

    """combine Force and Distance data (till 1st step and after last step)"""


    if TOMATO_param == 0:
        #find start of the first step and create a np array containing the data for fitting of the ds region (all data till first step)
        start_step1 = np.where(derivative_array[:, 1] == first_step_start)
        start_step1 = start_step1[0][0]

        raw_f_fitting_region_low = Force_Distance[0:int(len(F_low)/3), 0]
        raw_d_fitting_region_low = Force_Distance[0:int(len(F_low)/3), 1]

        #find end of the last step
        start_fitting_region = np.where(derivative_array[:, 1] == i_start)
        end_fitting_region = np.where(derivative_array[:, 1] == i_end)
        start_fitting_region = start_fitting_region[0][0]
        end_fitting_region = end_fitting_region[0][0]

        raw_f_fitting_region_up = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 0]
        raw_d_fitting_region_up = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 1]

    elif TOMATO_param == 1:
        start_step1 = find_nearest(Force_Distance[:, 1], i_start)
        raw_f_fitting_region_low = Force_Distance[0:(int(start_step1/3)), 0]
        raw_d_fitting_region_low = Force_Distance[0:(int(start_step1/3)), 1]

        start_fitting_region = find_nearest(Force_Distance[:, 1], i_start)
        end_fitting_region = find_nearest(Force_Distance[:, 1], i_end)
        raw_f_fitting_region_up = Force_Distance[start_fitting_region:end_fitting_region, 0]
        raw_d_fitting_region_up = Force_Distance[start_fitting_region:end_fitting_region, 1]

    d_fitting_region_ds, f_fitting_region_ds = raw_d_fitting_region_low, raw_f_fitting_region_low 

    Force_combined = np.concatenate((raw_f_fitting_region_low,raw_f_fitting_region_up))
    Distance_combined = np.concatenate((raw_d_fitting_region_low,raw_d_fitting_region_up)) 

    # downsample the data used for fitting to around 200 datapoints
    if len(Force_combined) > 200:
        Force_data =Force_combined[::int(len(Force_combined) / 200)]
        Distance_data = Distance_combined[::int(len(Force_combined) / 200)]

    else:
        Force_data = Force_combined
        Distance_data = Distance_combined




    """Fitting"""

    model_FU = lk.ewlc_odijk_distance("ds_part").subtract_independent_offset() + lk.ewlc_odijk_distance("RNA").subtract_independent_offset() + lk.distance_offset("offset") 
    model_FU = model_FU.invert()


    fit_FU=lk.FdFit(model_FU)

    fit_FU.add_data('Fully_unfolded', Force_data, Distance_data)

    # distance shift bounds
    fit_FU["offset/d_offset"].value = input_fitting['offset_d']
    fit_FU["offset/d_offset"].upper_bound = input_fitting['offset_d_up']
    fit_FU["offset/d_offset"].lower_bound = input_fitting['offset_d_low']
    fit_FU["offset/d_offset"].unit='nm'

    """ds handles part"""
    # Persistance length bounds
    fit_FU["ds_part/Lp"].value = input_fitting['lp_ds']
    fit_FU["ds_part/Lp"].upper_bound = input_fitting['lp_ds_up']
    fit_FU["ds_part/Lp"].lower_bound = input_fitting['lp_ds_low']
   
    # Force shift bounds
    fit_FU["ds_part/f_offset"].value = input_fitting['offset_f']
    fit_FU["ds_part/f_offset"].upper_bound = input_fitting['offset_f_up']
    fit_FU["ds_part/f_offset"].lower_bound = input_fitting['offset_f_low']


    # stiffnes
    fit_FU["ds_part/St"].value = input_fitting['ds_stiff']
    fit_FU["ds_part/St"].upper_bound = input_fitting['ds_stiff_up']
    fit_FU["ds_part/St"].lower_bound = input_fitting['ds_stiff_low']

    # contour length
    Lc_initial_guess = input_fitting['lc_ds']  # nm
    Lc_range = 5
    fit_FU["ds_part/Lc"].upper_bound = Lc_initial_guess + Lc_range
    fit_FU["ds_part/Lc"].lower_bound = Lc_initial_guess - Lc_range
    fit_FU["ds_part/Lc"].value = Lc_initial_guess
    fit_FU["ds_part/Lc"].unit = 'nm'
    fit_FU["ds_part/Lc"].fixed = 'True'

    """ssRNA part - full length"""

    # Persistance length bounds
    fit_FU["RNA/Lp"].value = input_fitting['lp_ss']
    fit_FU["RNA/Lp"].lower_bound = 0.8
    fit_FU["RNA/Lp"].upper_bound = 2
    fit_FU["RNA/Lp"].fixed = 'True'




    # stiffnes
    fit_FU["RNA/St"].value = input_fitting['ss_stiff']
    fit_FU["RNA/St"].lower_bound = input_fitting['ss_stiff_low']
    fit_FU["RNA/St"].upper_bound = input_fitting['ss_stiff_up']
    fit_FU["RNA/St"].fixed = 'True'

    # contour length
    fit_FU["RNA/Lc"].value = input_fitting['lc_ss_up']
    fit_FU["RNA/Lc"].lower_bound = 0
    fit_FU["RNA/Lc"].upper_bound = input_fitting['lc_ss_up']
    fit_FU["RNA/Lc"].fixed = 'True'
    fit_FU["RNA/Lc"].unit = 'nm'

    fit_FU["RNA/f_offset"].value = fit_FU["ds_part/f_offset"].value
    fit_FU["RNA/f_offset"].lower_bound =  input_fitting['offset_f_low']
    fit_FU["RNA/f_offset"].upper_bound = input_fitting['offset_f_up']





    fit_FU.fit()
    fit_qual = fit_FU.log_likelihood()
    print(fit_FU.params)


    """recreate the ds model for plotting"""
    model_ds = lk.ewlc_odijk_distance("ds_part").subtract_independent_offset() + lk.distance_offset("offset") 
    model_ds = model_ds.invert()

    fit_ds=lk.FdFit(model_ds)

    fit_ds.add_data('Fully_unfolded', Force_data, Distance_data)

    # distance shift bounds
    fit_ds["offset/d_offset"].value = fit_FU["offset/d_offset"].value
    # Persistance length bounds
    fit_ds["ds_part/Lp"].value = fit_FU["ds_part/Lp"].value
    # Force shift bounds
    fit_ds["ds_part/f_offset"].value = fit_FU["ds_part/f_offset"].value
    # stiffnes
    fit_ds["ds_part/St"].value = fit_FU["ds_part/St"].value
    # contour length
    fit_ds["ds_part/Lc"].value = fit_FU["ds_part/Lc"].value


    # calculate the integral until the first unfolding step
    # used to calculate the work done by the machine
    distance_integral = np.arange(min(Force_Distance[:, 1]), i_start)
    ds_integral = model_FU(distance_integral, fit_FU.params)
    area_ds = simps(ds_integral)
    print("area_ds = " + str(area_ds))

    # export the fitting parameters
    ds_fit_dict = {
        'filename': filename_i,
        'model': 'WLC',
        'model_ss': fit_FU,
        #'fit_model': fit_ds,
        'log_likelihood': fit_qual,
        'Lc_ds': fit_FU["ds_part/Lc"].value,
        'Lp_ds': fit_FU["ds_part/Lp"].value,
        'Lp_ds_stderr': fit_FU["ds_part/Lp"].stderr,
        'St_ds': fit_FU["ds_part/St"].value,
        'f_offset': fit_FU["ds_part/f_offset"].value,
        'd_offset': fit_FU["offset/d_offset"].value
    }
    print(ds_fit_dict)
    return ds_fit_dict, area_ds, start_step1




def fitting_FU_ss(filename_i, input_settings, export_data, input_fitting, i_start, i_end, Force_Distance, fix, max_range, derivative_array, F_low, TOMATO_param):
    global model_ss
    global ss_fit_dict
    #print("beginning of ss folding reached")
    if TOMATO_param == 0:
        end_fitting_region = np.where(derivative_array[:, 1] == i_end)
        end_fitting_region = end_fitting_region[0][0]
        raw_f_fitting_region = Force_Distance[0:end_fitting_region * input_settings['step_d'] + len(F_low), 0]
        raw_d_fitting_region = Force_Distance[0:end_fitting_region * input_settings['step_d'] + len(F_low), 1]
    elif TOMATO_param == 1:
        end_fitting_region = find_nearest(Force_Distance[:, 1], i_end)
        raw_f_fitting_region = Force_Distance[0:end_fitting_region, 0]
        raw_d_fitting_region = Force_Distance[0:end_fitting_region, 1]

    # downsample the data used for fitting to around 200 datapoints
    if len(raw_f_fitting_region) > 200:
        f_fitting_region_ss = raw_f_fitting_region[::int(len(raw_f_fitting_region) / 200)]
        d_fitting_region_ss = raw_d_fitting_region[::int(len(raw_f_fitting_region) / 200)]
    else:
        f_fitting_region_ss = raw_f_fitting_region
        d_fitting_region_ss = raw_d_fitting_region

    #if input_fitting['WLC+FJC'] == 1:
     #   model_ss = lk.odijk("DNA_2") + lk.freely_jointed_chain("RNA")
    #elif input_fitting['WLC+WLC'] == 1:
     #   model_ss = lk.odijk("DNA_2") + lk.odijk("RNA")
    
    if input_fitting['WLC+FJC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.efjc_distance("RNA")
    elif input_fitting['WLC+WLC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.ewlc_odijk_distance("RNA") 
    model_ss = model_ss.invert().subtract_independent_offset()+ lk.force_offset("DNA")
    fit_ss = lk.FdFit(model_ss)
    #print("ss model created")
    fit_ss.add_data("ss_part", f_fitting_region_ss, d_fitting_region_ss)
    #print("data added to ss model")
    # ds part parameters
    # Persistance length bounds

    # Lp_ds_range=fit_ds["DNA/Lp"].value/10
    fit_ss["DNA_2/Lp"].value = ds_fit_dict['Lp_ds']
    fit_ss["DNA_2/Lp"].upper_bound = ds_fit_dict['Lp_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lp"].lower_bound = ds_fit_dict['Lp_ds'] * (1 - max_range / 100)
    # if fix==1:
    fit_ss["DNA_2/Lp"].fixed = 'True'

    fit_ss["DNA/f_offset"].upper_bound = 5
    fit_ss["DNA/f_offset"].lower_bound = -5
    fit_ss["DNA/f_offset"].value = ds_fit_dict['f_offset']
    fit_ss["DNA/f_offset"].fixed = 'True'

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = ds_fit_dict['d_offset']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = 'True'

    # contour length
    # Lc_ds_range=Lc_initial_guess/100 # nm
    fit_ss["DNA_2/Lc"].upper_bound = ds_fit_dict['Lc_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = ds_fit_dict['Lc_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = ds_fit_dict['Lc_ds']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    # if fix==1:
    fit_ss["DNA_2/Lc"].fixed = 'True'

    # stifness

    fit_ss["DNA_2/St"].upper_bound = ds_fit_dict['St_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = ds_fit_dict['St_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = ds_fit_dict['St_ds']
    if fix == 1:
        fit_ss["DNA_2/St"].fixed = 'True'

    # ss part parameters

    # Persistance length bounds
    fit_ss["RNA/Lp"].value = input_fitting['lp_ss']
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1:
        fit_ss["RNA/Lp"].fixed = 'True'

    # stiffnes
    fit_ss["RNA/St"].value = input_fitting['ss_stiff']
    fit_ss["RNA/St"].lower_bound = input_fitting['ss_stiff_low']
    fit_ss["RNA/St"].upper_bound = input_fitting['ss_stiff_up']
    if fix == 1:
        fit_ss["RNA/St"].fixed = 'True'
    # contour length
    fit_ss["RNA/Lc"].value = input_fitting['lc_ss']
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].upper_bound = input_fitting['lc_ss_up']

    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()
    print(fit_ss.params)

    # calculate the integrals of the fitted functions

    distance_integral = np.arange(min(Force_Distance[:, 1]), i_end)
    ds_integral = model_ds(distance_integral, fit_ds.params)
    area_ds = simps(ds_integral)

    fit_qual = fit_ss.log_likelihood()

    if input_fitting["WLC+WLC"] == 1:
        fitting_model = "WLC+WLC"
    elif input_fitting["WLC+FJC"] == 1:
        fitting_model = "WLC+FJC"

    ss_fit_dict = {
        'filename': filename_i,
        'model': fitting_model,
        'model_ss': fit_ss,
        'log_likelihood': fit_qual,
        'Lc_ds': fit_ss["DNA_2/Lc"].value,
        'Lp_ds': fit_ss["DNA_2/Lp"].value,
        'St_ds': fit_ss["DNA_2/St"].value,
        'Lc_ss': fit_ss["RNA/Lc"].value,
        'Lc_ss_stderr': fit_ss["RNA/Lc"].stderr,
        'Lp_ss': fit_ss["RNA/Lp"].value,
        'St_ss': fit_ss["RNA/St"].value,
        'f_offset': fit_ss["DNA/f_offset"].value,
        'd_offset': fit_ss["inv(DNA_2_with_RNA)/d_offset"].value
    }

    return fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, area_ds










def fitting_ss(filename_i, input_settings, export_data, input_fitting, i_start, i_end, Force_Distance, fix, max_range, derivative_array, F_low, TOMATO_param):
    global model_ss
    global ss_fit_dict
    #print("beginning of ss folding reached")
    if TOMATO_param == 0:
        start_fitting_region = np.where(derivative_array[:, 1] == i_start)
        end_fitting_region = np.where(derivative_array[:, 1] == i_end)
        start_fitting_region = start_fitting_region[0][0]
        end_fitting_region = end_fitting_region[0][0]
        raw_f_fitting_region = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 0]
        raw_d_fitting_region = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 1]
    elif TOMATO_param == 1:
        start_fitting_region = find_nearest(Force_Distance[:, 1], i_start)
        end_fitting_region = find_nearest(Force_Distance[:, 1], i_end)
        raw_f_fitting_region = Force_Distance[start_fitting_region:end_fitting_region, 0]
        raw_d_fitting_region = Force_Distance[start_fitting_region:end_fitting_region, 1]

    # downsample the data used for fitting to around 200 datapoints
    if len(raw_f_fitting_region) > 200:
        f_fitting_region_ss = raw_f_fitting_region[::int(len(raw_f_fitting_region) / 200)]
        d_fitting_region_ss = raw_d_fitting_region[::int(len(raw_f_fitting_region) / 200)]
    else:
        f_fitting_region_ss = raw_f_fitting_region
        d_fitting_region_ss = raw_d_fitting_region

    #if input_fitting['WLC+FJC'] == 1:
     #   model_ss = lk.odijk("DNA_2") + lk.freely_jointed_chain("RNA")
    #elif input_fitting['WLC+WLC'] == 1:
     #   model_ss = lk.odijk("DNA_2") + lk.odijk("RNA")
    
    if input_fitting['WLC+FJC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.efjc_distance("RNA")
    elif input_fitting['WLC+WLC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.ewlc_odijk_distance("RNA") 
    model_ss = model_ss.invert().subtract_independent_offset()+ lk.force_offset("DNA")
    fit_ss = lk.FdFit(model_ss)
    #print("ss model created")
    fit_ss.add_data("ss_part", f_fitting_region_ss, d_fitting_region_ss)
    #print("data added to ss model")
    # ds part parameters
    # Persistance length bounds

    # Lp_ds_range=fit_ds["DNA/Lp"].value/10
    fit_ss["DNA_2/Lp"].value = ds_fit_dict['Lp_ds']
    fit_ss["DNA_2/Lp"].upper_bound = ds_fit_dict['Lp_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lp"].lower_bound = ds_fit_dict['Lp_ds'] * (1 - max_range / 100)
    # if fix==1:
    fit_ss["DNA_2/Lp"].fixed = 'True'

    fit_ss["DNA/f_offset"].upper_bound = 5
    fit_ss["DNA/f_offset"].lower_bound = -5
    fit_ss["DNA/f_offset"].value = ds_fit_dict['f_offset']
    fit_ss["DNA/f_offset"].fixed = 'True'

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = ds_fit_dict['d_offset']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = 'True'

    # contour length
    # Lc_ds_range=Lc_initial_guess/100 # nm
    fit_ss["DNA_2/Lc"].upper_bound = ds_fit_dict['Lc_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = ds_fit_dict['Lc_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = ds_fit_dict['Lc_ds']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    # if fix==1:
    fit_ss["DNA_2/Lc"].fixed = 'True'

    # stifness

    fit_ss["DNA_2/St"].upper_bound = ds_fit_dict['St_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = ds_fit_dict['St_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = ds_fit_dict['St_ds']
    if fix == 1:
        fit_ss["DNA_2/St"].fixed = 'True'

    # ss part parameters

    # Persistance length bounds
    fit_ss["RNA/Lp"].value = input_fitting['lp_ss']
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1:
        fit_ss["RNA/Lp"].fixed = 'True'

    # stiffnes
    fit_ss["RNA/St"].value = input_fitting['ss_stiff']
    fit_ss["RNA/St"].lower_bound = input_fitting['ss_stiff_low']
    fit_ss["RNA/St"].upper_bound = input_fitting['ss_stiff_up']
    if fix == 1:
        fit_ss["RNA/St"].fixed = 'True'
    # contour length
    fit_ss["RNA/Lc"].value = input_fitting['lc_ss']
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].upper_bound = input_fitting['lc_ss_up']

    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()
    print(fit_ss.params)

    # calculate the integrals of the fitted functions
    distance_integral_fit_start = np.arange(min(Force_Distance[:, 1]), i_start)
    ss_integral_start = model_ss(distance_integral_fit_start, fit_ss.params)
    area_ss_fit_start = simps(ss_integral_start)
    print("area_ss_start = " + str(area_ss_fit_start))

    distance_integral_fit_end = np.arange(min(Force_Distance[:, 1]), i_end)
    ss_integral_end = model_ss(distance_integral_fit_end, fit_ss.params)
    area_ss_fit_end = simps(ss_integral_end)
    print("area_ss_end = " + str(area_ss_fit_end))

    fit_qual = fit_ss.log_likelihood()

    if input_fitting["WLC+WLC"] == 1:
        fitting_model = "WLC+WLC"
    elif input_fitting["WLC+FJC"] == 1:
        fitting_model = "WLC+FJC"

    ss_fit_dict = {
        'filename': filename_i,
        'model': fitting_model,
        'model_ss': fit_ss,
        'model_ss_TOMATO': model_ss,
        'log_likelihood': fit_qual,
        'Lc_ds': fit_ss["DNA_2/Lc"].value,
        'Lp_ds': fit_ss["DNA_2/Lp"].value,
        'St_ds': fit_ss["DNA_2/St"].value,
        'Lc_ss': fit_ss["RNA/Lc"].value,
        'Lc_ss_stderr': fit_ss["RNA/Lc"].stderr,
        'Lp_ss': fit_ss["RNA/Lp"].value,
        'St_ss': fit_ss["RNA/St"].value,
        'f_offset': fit_ss["DNA/f_offset"].value,
        'd_offset': fit_ss["inv(DNA_2_with_RNA)/d_offset"].value
    }

    return fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, area_ss_fit_start, area_ss_fit_end

def plot_fit(fit, start_force_ss, start_distance_ss, Force_Distance, save_folder, filename_i, start_time):
    distance = np.arange(min(Force_Distance[:, 1]), max(Force_Distance[:, 1]) + 50, 2)
    F_ds_model = model_ds(distance, fit_ds.params)
    
    legend_elements = [
        Line2D([0], [0], color='k', lw=1, alpha=0.85),
        Line2D([0], [0], color='gray', linestyle='dashed', lw=1)
    ]

    diff_colors = ['b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm']

    font_size = 20
    line_thickness = 2
    min_x_value, max_x_value = min(Force_Distance[:, 1])-10, max(Force_Distance[:, 1])+10
    min_y_value, max_y_value = -1 , max(Force_Distance[:, 0])+3
    Plot_title=filename_i

    plt.plot(Force_Distance[:, 1], Force_Distance[:, 0], 'k', alpha=0.85)
    plt.axis([min(Force_Distance[:, 1]) - 50, max(Force_Distance[:, 1]) + 50, 0, max(Force_Distance[:, 0]) + 15])
    plt.scatter(d_fitting_region_ds, f_fitting_region_ds, color=diff_colors[0], s=4)
    plt.plot(distance, F_ds_model, linestyle='dashed', color=diff_colors[0], linewidth=0.5, alpha=0.85)
    plt.xlabel('Relative distance, nm', fontsize=font_size)
    plt.ylabel('Force, pN', fontsize=font_size)
    plt.legend(legend_elements, ['FD-Curve', 'Part used for fitting', 'Fitted WLC model'], fontsize=12)
    
    fit_data = {"distance": distance, "Fit_ds": F_ds_model}
    
    for i in range(len(fit)):
        F_ss_model = model_ss(distance, fit[i].params)
        plt.scatter(start_distance_ss[i], start_force_ss[i], s=4, color=diff_colors[i+1])
        plt.plot(distance, F_ss_model, linestyle='dashed', color=diff_colors[i+1], linewidth=0.5, alpha=0.85)
        fit_data[f"Fit_ss_{i+1}"] = F_ss_model
    

  


    plt.title(Plot_title, fontsize=12)

    # Set the tick and axis properties
    ax = plt.gca()  # Get current axis
    ax.spines['top'].set_linewidth(2)    # Set the width of the top spine
    ax.spines['right'].set_linewidth(2)  # Set the width of the right spine
    ax.spines['left'].set_linewidth(2)   # Set the width of the left spine
    ax.spines['bottom'].set_linewidth(2) # Set the width of the bottom spine
    ax.spines['top'].set_color('white') 
    ax.spines['right'].set_color('white')

    # Add a scale bar
    scalebar_length = 100  # Length of the scale bar in x-axis units
    scalebar_height = 0.01 * (max_y_value - min_y_value)  # Relative height of the scale bar
    scalebar_x_position = max_x_value-scalebar_length*1.001  # Position of the scale bar on the x-axis
    scalebar_y_position = min_y_value + scalebar_height * 2  # Position on the y-axis


    # Create a rectangle patch as a scale bar
    scalebar = patches.Rectangle((scalebar_x_position, scalebar_y_position), scalebar_length, scalebar_height, color='black')

    # Add the patch to the axes
    plt.gca().add_patch(scalebar)


    # Customizing tick parameters
    plt.tick_params(axis='both', which='major', labelsize=font_size, direction='in', length=6, width=2)
    #plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1)
    # Hide x-axis labels and ticks
    plt.tick_params(axis='x', which='both', labelbottom=False, direction='in', length=6, width=2)
    # Optionally add text to label the scale bar
    plt.text(scalebar_x_position+scalebar_length/2, scalebar_y_position + scalebar_height * 2, f'{scalebar_length} nm', fontsize=font_size, ha='center')


    # Set axis limits
    plt.xlim(min_x_value, max_x_value)  # Set the range for the x-axis
    plt.ylim(min_y_value, max_y_value)  # Set the range for the y-axis



    # Save plot
    plotname = f"{save_folder}/{filename_i}_fit_{start_time}.png"
    plotname_svg = f"{save_folder}/{filename_i}_fit_{start_time}.svg"
    plt.savefig(plotname, dpi=600)
    plt.savefig(plotname_svg, format='svg')
    plt.clf()
    
    # Save fit data to CSV
    fit_df = pd.DataFrame(fit_data)
    csv_filename = f"{save_folder}/{filename_i}_fit_data_{start_time}.csv"
    fit_df.to_csv(csv_filename, index=False)
    print(f"Fit data saved to {csv_filename}")

"""
def plot_fit(fit, start_force_ss, start_distance_ss, Force_Distance, save_folder, filename_i, start_time):
    distance = np.arange(min(Force_Distance[:, 1]), max(Force_Distance[:, 1]) + 50, 2)
    F_ds_model = model_ds(distance, fit_ds.params)

    legend_elements = [
        Line2D([0], [0], color='k', lw=1, alpha=0.85),
        # Line2D([0], [0], color='r', lw=1),
        Line2D([0], [0], color='gray', linestyle='dashed', lw=1)
    ]

    diff_colors = ['b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm']

    plt.plot(Force_Distance[:, 1], Force_Distance[:, 0], 'k', alpha=0.85)
    plt.axis([min(Force_Distance[:, 1]) - 50, max(Force_Distance[:, 1]) + 50, 0, max(Force_Distance[:, 0]) + 15])
    plt.scatter(d_fitting_region_ds, f_fitting_region_ds, color=diff_colors[0], s=4)
    plt.plot(distance, F_ds_model, linestyle='dashed', color=diff_colors[0], linewidth=0.5, alpha=0.85)
    plt.ylabel("Force [pN]")
    plt.xlabel("Distance [nm]")
    plt.legend(legend_elements, ['FD-Curve', 'Part used for fitting', 'Fitted WLC model'])
    
    for i in range(0, len(fit)):
        F_ss_model = model_ss(distance, fit[i].params)
        plt.scatter(start_distance_ss[i], start_force_ss[i], s=4, color=diff_colors[i+1])
        plt.plot(distance, F_ss_model, linestyle='dashed', color=diff_colors[i+1], linewidth=0.5, alpha=0.85)

    plotname = save_folder + "/" + filename_i + "_fit_" + start_time + ".png"

    plt.savefig(plotname, dpi=600)
    plt.clf()
"""