"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


import matplotlib.pyplot as plt
import lumicks.pylake as lk
import numpy as np
import pandas as pd
from scipy.integrate import simps
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import pickle
import os
import traceback


"""define the functions used for fitting"""


# Helper function to safely get parameter values from FdFit objects (supports inverted models)
def safe_get_param(fit_obj, param_name):
    """
    Attempts to retrieve a parameter value and stderr from a fit object.
    Handles fixed parameters (where stderr is None) and inverted model key differences.
    """
    # Try direct lookup first
    try:
        if param_name in fit_obj.params:
            p = fit_obj.params[param_name]
            val = p.value if (hasattr(p, 'value') and p.value is not None) else np.nan
            err = p.stderr if (hasattr(p, 'stderr') and p.stderr is not None) else np.nan
            unit = getattr(p, 'unit', '') or ''
            return val, err, unit
    except Exception:
        pass
    
    # Fallback: iterate and match by suffix (for inverted model wrappers)
    try:
        suffix = param_name.split('/')[-1]
        for key, p in fit_obj.params.items():
            if key == param_name or key.endswith('/' + suffix):
                val = p.value if (hasattr(p, 'value') and p.value is not None) else np.nan
                err = p.stderr if (hasattr(p, 'stderr') and p.stderr is not None) else np.nan
                unit = getattr(p, 'unit', '') or ''
                return val, err, unit
    except Exception:
        pass

    return np.nan, np.nan, ''

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
    #print(fit_ds.params)

    # calculate the integral until the first unfolding step
    distance_integral = np.arange(min(Force_Distance[:, 1]), i_start)
    ds_integral = model_ds(distance_integral, fit_ds.params)
    area_ds = simps(ds_integral)
    #print("area_ds = " + str(area_ds))

    # Export parameters safely
    def get_val(key):
        val, err, unit = safe_get_param(fit_ds, key)
        return val if val is not None else 0, err if err is not None else np.nan

    ds_fit_dict = {
        'log_likelihood': fit_qual,
        'Lc_ds': get_val('ds_part/Lc')[0],
        'Lc_ds_stderr': get_val('ds_part/Lc')[1],
        'Lp_ds': get_val('ds_part/Lp')[0],
        'Lp_ds_stderr': get_val('ds_part/Lp')[1],
        'St_ds': get_val('ds_part/St')[0],
        'St_ds_stderr': get_val('ds_part/St')[1],
        'f_offset_ds': get_val('ds_part/f_offset')[0],
        'f_offset_ds_stderr': get_val('ds_part/f_offset')[1],
        'd_offset_ds': get_val('ds_part/d_offset')[0],
        'd_offset_ds_stderr': get_val('ds_part/d_offset')[1],
        'Lc_ss': 0,
        'Lc_ss_stderr': np.nan,
        'Lp_ss': np.nan,
        'Lp_ss_stderr': np.nan,
        'St_ss': np.nan,
        'St_ss_stderr': np.nan,
        'f_offset_ss': np.nan,
        'f_offset_ss_stderr': np.nan,
        'd_offset_ss': np.nan,
        'd_offset_ss_stderr': np.nan,
        'model_type': 'WLC_ODIJK_DS',
        'fit_status': 'success'
    }
    #print(ds_fit_dict)
    return ds_fit_dict, area_ds, start_step1, fit_ds

def fitting_FU(filename_i, input_settings, export_data, input_fitting,
               first_step_start, i_start, i_end, Force_Distance,
               derivative_array, F_low, TOMATO_param):
    global model_ss, ss_fit_dict
    global model_ds, fit_ds
    global ds_fit_dict
    global f_fitting_region_ds, d_fitting_region_ds
    global export_fit_ds
    global fitting_model

    # ----- Build fitting regions (unchanged logic) -----
    if TOMATO_param == 0:
        start_step1 = np.where(derivative_array[:, 1] == first_step_start)[0][0]

        raw_f_fitting_region_low = Force_Distance[0:int(len(F_low)/3), 0]
        raw_d_fitting_region_low = Force_Distance[0:int(len(F_low)/3), 1]

        start_fitting_region = np.where(derivative_array[:, 1] == i_start)[0][0]
        end_fitting_region   = np.where(derivative_array[:, 1] == i_end)[0][0]

        s = start_fitting_region * input_settings['step_d'] + len(F_low)
        e = end_fitting_region   * input_settings['step_d'] + len(F_low)

        raw_f_fitting_region_up = Force_Distance[s:e, 0]
        raw_d_fitting_region_up = Force_Distance[s:e, 1]
        i_start_distance = Force_Distance[s, 1] if s < len(Force_Distance) else Force_Distance[-1, 1]

    else:  # TOMATO_param == 1
        start_step1 = find_nearest(Force_Distance[:, 1], i_start)

        raw_f_fitting_region_low = Force_Distance[0:(int(start_step1/3)), 0]
        raw_d_fitting_region_low = Force_Distance[0:(int(start_step1/3)), 1]

        s = find_nearest(Force_Distance[:, 1], i_start)
        e = find_nearest(Force_Distance[:, 1], i_end)

        raw_f_fitting_region_up = Force_Distance[s:e, 0]
        raw_d_fitting_region_up = Force_Distance[s:e, 1]
        i_start_distance = Force_Distance[s, 1]

    d_fitting_region_ds, f_fitting_region_ds = raw_d_fitting_region_low, raw_f_fitting_region_low

    Force_combined    = np.concatenate((raw_f_fitting_region_low, raw_f_fitting_region_up))
    Distance_combined = np.concatenate((raw_d_fitting_region_low, raw_d_fitting_region_up))

    # ----- Downsample safely to ~200 points -----
    if len(Force_combined) > 200:
        step = max(int(np.ceil(len(Force_combined) / 200)), 1)
        Force_data    = Force_combined[::step]
        Distance_data = Distance_combined[::step]
    else:
        Force_data    = Force_combined
        Distance_data = Distance_combined

    # ----- Clean data -----
    m = np.isfinite(Force_data) & np.isfinite(Distance_data)
    Force_data, Distance_data = Force_data[m], Distance_data[m]
    if Force_data.size == 0:
        raise ValueError("No finite data in fitting region.")

    # ----- Build inverted model: distance -> force -----
    model_FU = (
        lk.ewlc_odijk_distance("ds_part").subtract_independent_offset() +
        lk.ewlc_odijk_distance("RNA").subtract_independent_offset() +
        lk.distance_offset("offset")
    ).invert()

    fit_FU = lk.FdFit(model_FU)

    fit_FU.add_data('Fully_unfolded', Force_data, Distance_data)

    # ----- Parameter initialization -----
    # Distance offset
    fit_FU["offset/d_offset"].value       = input_fitting['offset_d']
    fit_FU["offset/d_offset"].upper_bound = input_fitting['offset_d_up']
    fit_FU["offset/d_offset"].lower_bound = input_fitting['offset_d_low']
    fit_FU["offset/d_offset"].unit        = 'nm'

    fit_FU["ds_part/Lp"].value       = max(float(input_fitting['lp_ds']),      1e-6)
    fit_FU["ds_part/Lp"].upper_bound = max(float(input_fitting['lp_ds_up']),   1e-6)
    fit_FU["ds_part/Lp"].lower_bound = max(float(input_fitting['lp_ds_low']),  1e-6)

    fit_FU["ds_part/St"].value       = max(float(input_fitting['ds_stiff']),     1e-6)
    fit_FU["ds_part/St"].upper_bound = max(float(input_fitting['ds_stiff_up']),  1e-6)
    fit_FU["ds_part/St"].lower_bound = max(float(input_fitting['ds_stiff_low']), 1e-6)

    Lc0 = float(input_fitting['lc_ds'])
    Lc_rng = 5.0
    fit_FU["ds_part/Lc"].value       = max(Lc0, 1e-6)
    fit_FU["ds_part/Lc"].upper_bound = max(Lc0 + Lc_rng, 1e-6)
    fit_FU["ds_part/Lc"].lower_bound = max(Lc0 - Lc_rng, 1e-6)
    fit_FU["ds_part/Lc"].unit        = 'nm'
    fit_FU["ds_part/Lc"].fixed       = True

    fit_FU["RNA/Lp"].value       = float(input_fitting['lp_ss'])
    fit_FU["RNA/Lp"].lower_bound = 0.8
    fit_FU["RNA/Lp"].upper_bound = 2.0
    fit_FU["RNA/Lp"].fixed       = True

    fit_FU["RNA/St"].value       = float(input_fitting['ss_stiff'])
    fit_FU["RNA/St"].lower_bound = float(input_fitting['ss_stiff_low'])
    fit_FU["RNA/St"].upper_bound = float(input_fitting['ss_stiff_up'])
    fit_FU["RNA/St"].fixed       = True

    fit_FU["RNA/Lc"].value       = float(input_fitting['lc_ss_up'])
    fit_FU["RNA/Lc"].lower_bound = 0.0
    fit_FU["RNA/Lc"].upper_bound = float(input_fitting['lc_ss_up'])
    fit_FU["RNA/Lc"].unit        = 'nm'
    fit_FU["RNA/Lc"].fixed       = True

    # Force offset handling
    fmin = float(np.min(Force_data))
    eps = 1e-6
    start_f_offset = float(input_fitting['offset_f'])
    safe_upper = min(float(input_fitting['offset_f_up']), fmin - eps)
    safe_lower = float(input_fitting['offset_f_low'])

    if safe_upper <= safe_lower:
        safe_upper = fmin - eps
        safe_lower = safe_upper - max(0.1 * max(fmin, 1.0), 1e-3)

    start_f_offset = min(start_f_offset, safe_upper - 0.1 * abs(safe_upper - safe_lower))

    fit_FU["ds_part/f_offset"].lower_bound = safe_lower
    fit_FU["ds_part/f_offset"].upper_bound = safe_upper
    fit_FU["ds_part/f_offset"].value       = start_f_offset

    fit_FU["RNA/f_offset"].lower_bound = safe_lower
    fit_FU["RNA/f_offset"].upper_bound = safe_upper
    fit_FU["RNA/f_offset"].value       = start_f_offset

    if (fmin - start_f_offset) <= 0:
        fit_FU["ds_part/f_offset"].value -= max(0.1 * max(fmin, 1.0), 1e-3)
        fit_FU["RNA/f_offset"].value     = fit_FU["ds_part/f_offset"].value

    fit_FU.fit()
    fit_qual = fit_FU.log_likelihood()
    #print(fit_FU.params)

    # ----- ds-only model for plotting -----
    model_ds = (lk.ewlc_odijk_distance("ds_part").subtract_independent_offset()
                + lk.distance_offset("offset")).invert()
    fit_ds = lk.FdFit(model_ds)
    fit_ds.add_data('Fully_unfolded', Distance_data, Force_data)

    fit_ds["offset/d_offset"].value   = fit_FU["offset/d_offset"].value
    fit_ds["ds_part/Lp"].value        = fit_FU["ds_part/Lp"].value
    fit_ds["ds_part/f_offset"].value  = fit_FU["ds_part/f_offset"].value
    fit_ds["ds_part/St"].value        = fit_FU["ds_part/St"].value
    fit_ds["ds_part/Lc"].value        = fit_FU["ds_part/Lc"].value

    d_min = float(np.nanmin(Force_Distance[:, 1]))
    distance_integral = np.linspace(d_min, float(i_start_distance), 400)
    ds_integral = model_FU(distance_integral, fit_FU.params)
    area_ds = simps(ds_integral, distance_integral)
    #print("area_ds =", area_ds)

    # Safe extraction for inverted model keys
    def get_fit_val(key):
        val, err, _ = safe_get_param(fit_FU, key)
        return val if val is not None else 0, err if err is not None else np.nan

    ds_fit_dict = {
        'model': 'WLC',
        'log_likelihood': fit_qual,
        'Lc_ds': get_fit_val('ds_part/Lc')[0],
        'Lc_ds_stderr': get_fit_val('ds_part/Lc')[1],
        'Lp_ds': get_fit_val('ds_part/Lp')[0],
        'Lp_ds_stderr': get_fit_val('ds_part/Lp')[1],
        'St_ds': get_fit_val('ds_part/St')[0],
        'St_ds_stderr': get_fit_val('ds_part/St')[1],
        'f_offset_ds': get_fit_val('ds_part/f_offset')[0],
        'f_offset_ds_stderr': get_fit_val('ds_part/f_offset')[1],
        'd_offset_ds': get_fit_val('offset/d_offset')[0],
        'd_offset_ds_stderr': get_fit_val('offset/d_offset')[1],
        'Lc_ss_RNA': get_fit_val('RNA/Lc')[0],
        'Lc_ss_RNA_stderr': get_fit_val('RNA/Lc')[1],
        'Lp_ss_RNA': get_fit_val('RNA/Lp')[0],
        'Lp_ss_RNA_stderr': get_fit_val('RNA/Lp')[1],
        'St_ss_RNA': get_fit_val('RNA/St')[0],
        'St_ss_RNA_stderr': get_fit_val('RNA/St')[1],
        'f_offset_RNA': get_fit_val('RNA/f_offset')[0],
        'f_offset_RNA_stderr': get_fit_val('RNA/f_offset')[1],
        'model_type': 'WLC_ODIJK_FULLY_UNFOLDED',
        'fit_status': 'success',
        'Lc_ss': 0, 'Lc_ss_stderr': np.nan, 'Lp_ss': np.nan, 'Lp_ss_stderr': np.nan,
        'St_ss': np.nan, 'St_ss_stderr': np.nan, 'f_offset_ss': np.nan, 'f_offset_ss_stderr': np.nan,
        'd_offset_ss': np.nan, 'd_offset_ss_stderr': np.nan
    }
    return ds_fit_dict, area_ds, start_step1, fit_FU


def fitting_FU_ss(filename_i, input_settings, export_data, input_fitting, i_start, i_end, Force_Distance, fix, max_range, derivative_array, F_low, TOMATO_param):
    global model_ss
    global ss_fit_dict
    
    if TOMATO_param == 0:
        end_fitting_region = np.where(derivative_array[:, 1] == i_end)
        end_fitting_region = end_fitting_region[0][0]
        raw_f_fitting_region = Force_Distance[0:end_fitting_region * input_settings['step_d'] + len(F_low), 0]
        raw_d_fitting_region = Force_Distance[0:end_fitting_region * input_settings['step_d'] + len(F_low), 1]
    elif TOMATO_param == 1:
        end_fitting_region = find_nearest(Force_Distance[:, 1], i_end)
        raw_f_fitting_region = Force_Distance[0:end_fitting_region, 0]
        raw_d_fitting_region = Force_Distance[0:end_fitting_region, 1]

    if len(raw_f_fitting_region) > 200:
        f_fitting_region_ss = raw_f_fitting_region[::int(len(raw_f_fitting_region) / 200)]
        d_fitting_region_ss = raw_d_fitting_region[::int(len(raw_f_fitting_region) / 200)]
    else:
        f_fitting_region_ss = raw_f_fitting_region
        d_fitting_region_ss = raw_d_fitting_region

    if input_fitting['WLC+FJC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.efjc_distance("RNA")
    elif input_fitting['WLC+WLC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.ewlc_odijk_distance("RNA") 
    
    # CRITICAL FIX: Ensure we are using the correct inverted chain
    model_ss = model_ss.invert().subtract_independent_offset()+ lk.force_offset("DNA")
    
    fit_ss = lk.FdFit(model_ss)
    fit_ss.add_data("ss_part", f_fitting_region_ss, d_fitting_region_ss)

    # Parameters setup (unchanged logic)
    fit_ss["DNA_2/Lp"].value = ds_fit_dict['Lp_ds']
    fit_ss["DNA_2/Lp"].upper_bound = ds_fit_dict['Lp_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lp"].lower_bound = ds_fit_dict['Lp_ds'] * (1 - max_range / 100)
    if fix==1: fit_ss["DNA_2/Lp"].fixed = True

    fit_ss["DNA/f_offset"].upper_bound = 5
    fit_ss["DNA/f_offset"].lower_bound = -5
    fit_ss["DNA/f_offset"].value = ds_fit_dict['f_offset_ds']
    fit_ss["DNA/f_offset"].fixed = True

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = ds_fit_dict['d_offset_ds']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = True

    fit_ss["DNA_2/Lc"].upper_bound = ds_fit_dict['Lc_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = ds_fit_dict['Lc_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = ds_fit_dict['Lc_ds']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    if fix==1: fit_ss["DNA_2/Lc"].fixed = True

    fit_ss["DNA_2/St"].upper_bound = ds_fit_dict['St_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = ds_fit_dict['St_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = ds_fit_dict['St_ds']
    if fix == 1: fit_ss["DNA_2/St"].fixed = True

    fit_ss["RNA/Lp"].value = input_fitting['lp_ss']
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1: fit_ss["RNA/Lp"].fixed = True

    fit_ss["RNA/St"].value = input_fitting['ss_stiff']
    fit_ss["RNA/St"].lower_bound = input_fitting['ss_stiff_low']
    fit_ss["RNA/St"].upper_bound = input_fitting['ss_stiff_up']
    if fix == 1: fit_ss["RNA/St"].fixed = True

    fit_ss["RNA/Lc"].value = input_fitting['lc_ss']
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].upper_bound = input_fitting['lc_ss_up']
    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()
    #print(fit_ss.params)

    distance_integral = np.arange(min(Force_Distance[:, 1]), i_end)
    ds_integral = model_ds(distance_integral, fit_ds.params)
    area_ds = simps(ds_integral)

    fit_qual = fit_ss.log_likelihood()

    fitting_model = "WLC+WLC" if input_fitting["WLC+WLC"] == 1 else "WLC+FJC"

    # CRITICAL FIX: Use safe_get_param for inverted model keys
    def get_val(key):
        val, err, _ = safe_get_param(fit_ss, key)
        return val if val is not None else 0, err if err is not None else np.nan

    ss_fit_dict = {
        'log_likelihood': fit_qual,
        'Lc_ds': get_val('DNA_2/Lc')[0],
        'Lc_ds_stderr': get_val('DNA_2/Lc')[1],
        'Lp_ds': get_val('DNA_2/Lp')[0],
        'Lp_ds_stderr': get_val('DNA_2/Lp')[1],
        'St_ds': get_val('DNA_2/St')[0],
        'St_ds_stderr': get_val('DNA_2/St')[1],
        'f_offset_ds': get_val('DNA/f_offset')[0],
        'f_offset_ds_stderr': get_val('DNA/f_offset')[1],
        'd_offset_ds': get_val('inv(DNA_2_with_RNA)/d_offset')[0],
        'd_offset_ds_stderr': get_val('inv(DNA_2_with_RNA)/d_offset')[1],
        'Lc_ss': get_val('RNA/Lc')[0],
        'Lc_ss_stderr': get_val('RNA/Lc')[1],
        'Lp_ss': get_val('RNA/Lp')[0],
        'Lp_ss_stderr': get_val('RNA/Lp')[1],
        'St_ss': get_val('RNA/St')[0],
        'St_ss_stderr': get_val('RNA/St')[1],
        'f_offset_ss': get_val('DNA/f_offset')[0],
        'f_offset_ss_stderr': get_val('DNA/f_offset')[1],
        'd_offset_ss': get_val('inv(DNA_2_with_RNA)/d_offset')[0],
        'd_offset_ss_stderr': get_val('inv(DNA_2_with_RNA)/d_offset')[1],
        'model_type': fitting_model,
        'fit_status': 'success'
    }

    return fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, area_ds

def fitting_ss(filename_i, input_settings, export_data, input_fitting, i_start, i_end, Force_Distance, fix, max_range, derivative_array, F_low, TOMATO_param):
    global model_ss
    global ss_fit_dict
    
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

    if len(raw_f_fitting_region) > 200:
        f_fitting_region_ss = raw_f_fitting_region[::int(len(raw_f_fitting_region) / 200)]
        d_fitting_region_ss = raw_d_fitting_region[::int(len(raw_f_fitting_region) / 200)]
    else:
        f_fitting_region_ss = raw_f_fitting_region
        d_fitting_region_ss = raw_d_fitting_region

    if input_fitting['WLC+FJC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.efjc_distance("RNA")
    elif input_fitting['WLC+WLC'] == 1:
        model_ss = lk.ewlc_odijk_distance("DNA_2") + lk.ewlc_odijk_distance("RNA") 
    
    # CRITICAL FIX: Ensure inversion is applied correctly
    model_ss = model_ss.invert().subtract_independent_offset()+ lk.force_offset("DNA")
    
    fit_ss = lk.FdFit(model_ss)
    fit_ss.add_data("ss_part", f_fitting_region_ss, d_fitting_region_ss)

    # Parameters setup
    fit_ss["DNA_2/Lp"].value = ds_fit_dict['Lp_ds']
    fit_ss["DNA_2/Lp"].upper_bound = ds_fit_dict['Lp_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lp"].lower_bound = ds_fit_dict['Lp_ds'] * (1 - max_range / 100)
    if fix==1: fit_ss["DNA_2/Lp"].fixed = True

    fit_ss["DNA/f_offset"].upper_bound = 5
    fit_ss["DNA/f_offset"].lower_bound = -5
    fit_ss["DNA/f_offset"].value = ds_fit_dict['f_offset_ds']
    fit_ss["DNA/f_offset"].fixed = True

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = ds_fit_dict['d_offset_ds']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = True

    fit_ss["DNA_2/Lc"].upper_bound = ds_fit_dict['Lc_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = ds_fit_dict['Lc_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = ds_fit_dict['Lc_ds']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    if fix==1: fit_ss["DNA_2/Lc"].fixed = True

    fit_ss["DNA_2/St"].upper_bound = ds_fit_dict['St_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = ds_fit_dict['St_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = ds_fit_dict['St_ds']
    if fix == 1: fit_ss["DNA_2/St"].fixed = True

    fit_ss["RNA/Lp"].value = input_fitting['lp_ss']
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1: fit_ss["RNA/Lp"].fixed = True

    fit_ss["RNA/St"].value = input_fitting['ss_stiff']
    fit_ss["RNA/St"].lower_bound = input_fitting['ss_stiff_low']
    fit_ss["RNA/St"].upper_bound = input_fitting['ss_stiff_up']
    if fix == 1: fit_ss["RNA/St"].fixed = True

    fit_ss["RNA/Lc"].value = input_fitting['lc_ss']
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].upper_bound = input_fitting['lc_ss_up']
    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()
    #print(fit_ss.params)

    distance_integral_fit_start = np.arange(min(Force_Distance[:, 1]), i_start)
    ss_integral_start = model_ss(distance_integral_fit_start, fit_ss.params)
    area_ss_fit_start = simps(ss_integral_start)
    #print("area_ss_start = " + str(area_ss_fit_start))

    distance_integral_fit_end = np.arange(min(Force_Distance[:, 1]), i_end)
    ss_integral_end = model_ss(distance_integral_fit_end, fit_ss.params)
    area_ss_fit_end = simps(ss_integral_end)
    #print("area_ss_end = " + str(area_ss_fit_end))

    fit_qual = fit_ss.log_likelihood()
    fitting_model = "WLC+WLC" if input_fitting["WLC+WLC"] == 1 else "WLC+FJC"

    # CRITICAL FIX: Use safe_get_param for inverted model keys
    def get_val(key):
        val, err, _ = safe_get_param(fit_ss, key)
        return val if val is not None else 0, err if err is not None else np.nan

    ss_fit_dict = {
        'log_likelihood': fit_qual,
        'Lc_ds': get_val('DNA_2/Lc')[0],
        'Lc_ds_stderr': get_val('DNA_2/Lc')[1],
        'Lp_ds': get_val('DNA_2/Lp')[0],
        'Lp_ds_stderr': get_val('DNA_2/Lp')[1],
        'St_ds': get_val('DNA_2/St')[0],
        'St_ds_stderr': get_val('DNA_2/St')[1],
        'f_offset_ds': get_val('DNA/f_offset')[0],
        'f_offset_ds_stderr': get_val('DNA/f_offset')[1],
        'd_offset_ds': get_val('inv(DNA_2_with_RNA)/d_offset')[0],
        'd_offset_ds_stderr': get_val('inv(DNA_2_with_RNA)/d_offset')[1],
        'Lc_ss': get_val('RNA/Lc')[0],
        'Lc_ss_stderr': get_val('RNA/Lc')[1],
        'Lp_ss': get_val('RNA/Lp')[0],
        'Lp_ss_stderr': get_val('RNA/Lp')[1],
        'St_ss': get_val('RNA/St')[0],
        'St_ss_stderr': get_val('RNA/St')[1],
        'f_offset_ss': get_val('DNA/f_offset')[0],
        'f_offset_ss_stderr': get_val('DNA/f_offset')[1],
        'd_offset_ss': get_val('inv(DNA_2_with_RNA)/d_offset')[0],
        'd_offset_ss_stderr': get_val('inv(DNA_2_with_RNA)/d_offset')[1],
        'model_type': fitting_model,
        'fit_status': 'success'
    }
    return fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, area_ss_fit_start, area_ss_fit_end

def plot_fit(fit, start_force_ss, start_distance_ss, Force_Distance, Force_Distance_ds, save_folder, filename_i, start_time, export_data, model_FU=None, model_ds_final=None):
    """
    Updated plot_fit to also save the pyLake FdFit objects as .pkl files.
    """
    distance = np.arange(min(Force_Distance[:, 1]), max(Force_Distance[:, 1]) + 50, 2)
    
    # Ensure global variables are available for plotting if they exist
    try:
        F_ds_model = model_ds(distance, fit_ds.params) if 'model_ds' in globals() and 'fit_ds' in globals() else None
    except:
        F_ds_model = None

    legend_elements = [
        Line2D([0], [0], color='k', lw=1, alpha=0.5),
        Line2D([0], [0], color='gray', linestyle='dashed', lw=1)
    ]

    diff_colors = ['b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm']

    font_size = 20
    line_thickness = 2
    min_x_value, max_x_value = min(Force_Distance[:, 1])-10, max(Force_Distance[:, 1])+10
    min_y_value, max_y_value = -1 , max(Force_Distance[:, 0])+3
    Plot_title=filename_i

    plt.plot(Force_Distance[:, 1], Force_Distance[:, 0], 'k', alpha=1)
    plt.plot(Force_Distance_ds[:, 1], Force_Distance_ds[:, 0], 'k', alpha=0.2)
    plt.axis([min(Force_Distance[:, 1]) - 50, max(Force_Distance[:, 1]) + 50, 0, max(Force_Distance[:, 0]) + 15])
    
    try:
        plt.scatter(d_fitting_region_ds, f_fitting_region_ds, color=diff_colors[0],marker=".", s=6,  alpha=0.5)
        if F_ds_model is not None:
            plt.plot(distance, F_ds_model, linestyle='dashed', color=diff_colors[0], linewidth=1, alpha=1)
    except:
        pass

    plt.xlabel('Relative distance, nm', fontsize=font_size)
    plt.ylabel('Force, pN', fontsize=font_size)
    plt.legend(legend_elements, ['FD-Curve', 'Part used for fitting', 'Fitted WLC model'], fontsize=12)
    
    fit_data = {"distance": distance}
    if F_ds_model is not None:
        fit_data["Fit_ds"] = F_ds_model
    
    for i in range(len(fit)):
        try:
            F_ss_model = model_ss(distance, fit[i].params)
            plt.scatter(start_distance_ss[i], start_force_ss[i], s=6, marker=".",color=diff_colors[i+1], alpha=0.5)
            plt.plot(distance, F_ss_model, linestyle='dashed', color=diff_colors[i+1], linewidth=1, alpha=1)
            fit_data[f"Fit_ss_{i+1}"] = F_ss_model
        except Exception as e:
            print(f"Warning: Could not plot fit {i}: {e}")

    plt.title(Plot_title, fontsize=12)

    ax = plt.gca()
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_color('white') 
    ax.spines['right'].set_color('white')

    scalebar_length = 100
    scalebar_height = 0.01 * (max_y_value - min_y_value)
    scalebar_x_position = max_x_value-scalebar_length*1.001
    scalebar_y_position = min_y_value + scalebar_height * 2
    scalebar = patches.Rectangle((scalebar_x_position, scalebar_y_position), scalebar_length, scalebar_height, color='black')
    plt.gca().add_patch(scalebar)

    plt.tick_params(axis='both', which='major', labelsize=font_size, direction='in', length=6, width=2)
    plt.tick_params(axis='x', which='both', labelbottom=False, direction='in', length=6, width=2)
    plt.text(scalebar_x_position+scalebar_length/2, scalebar_y_position + scalebar_height * 2, f'{scalebar_length} nm', fontsize=font_size, ha='center')
    plt.xlim(min_x_value, max_x_value)
    plt.ylim(min_y_value, max_y_value)

    plotname = f"{save_folder}/{filename_i}_fit_{start_time}.png"
    plt.savefig(plotname, dpi=150)

    # Save SVG only if enabled
    if export_data.get('export_svg', True):  # Default to True
        plotname_svg = f"{save_folder}/{filename_i}_fit_{start_time}.svg"
        plt.savefig(plotname_svg, format='svg')
    #plt.clf()
    plt.close()
    
    fit_df = pd.DataFrame(fit_data)
    csv_filename = f"{save_folder}/{filename_i}_fit_data_{start_time}.csv"
    fit_df.to_csv(csv_filename, index=False)
    print(f"Fit data saved to {csv_filename}")

    # Save Pylake Models
    models_to_save = {}
    if model_ds_final is not None:
        models_to_save["ds_initial_fit"] = model_ds_final
    if model_FU is not None:
        models_to_save["full_unfolded_fit"] = model_FU
    if fit:
        for idx, fit_obj in enumerate(fit):
            models_to_save[f"ss_fit_{idx}"] = fit_obj
            
    if models_to_save:
        pkl_filename = f"{save_folder}/{filename_i}_pylake_models_{start_time}.pkl"
        try:
            with open(pkl_filename, 'wb') as f_out:
                pickle.dump(models_to_save, f_out, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"PyLake models saved to {pkl_filename}")
        except Exception as e:
            print(f"Error saving pickle file for {filename_i}: {e}")
    else:
        print("No PyLake models found to save.")