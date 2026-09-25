"""Copyright 2024 Lukáš Pekárek & Stefan Buck"""


import pandas as pd
import h5py
import numpy as np
from pathlib import Path
import lumicks.pylake as lk
import traceback
import gc

# relative imports
from FRIED_POTATO_fitting import fitting_ds, fitting_ss, plot_fit, fitting_FU, fitting_FU_ss
from FRIED_POTATO_preprocessing import preprocess_RAW, trim_data, create_derivative
from FRIED_POTATO_find_steps import find_steps_F,find_steps_F_old, find_steps_PD, find_steps_PD_old, find_common_steps, calc_integral, save_figure
from FRIED_POTATO_processMultiH5 import split_H5


# ---------------------------------------------------------------------------
# Canonical fit-result schema. MUST match the header written to
# total_results_<timestamp>.csv in start_subprocess().
# ---------------------------------------------------------------------------
FIT_HEADER = [
    'model_type', 'log_likelihood',
    'Lc_ds', 'Lc_ds_stderr', 'Lp_ds', 'Lp_ds_stderr',
    'St_ds', 'St_ds_stderr', 'f_offset_ds', 'f_offset_ds_stderr',
    'd_offset_ds', 'd_offset_ds_stderr',
    'Lc_ss', 'Lc_ss_stderr', 'Lp_ss', 'Lp_ss_stderr',
    'St_ss', 'St_ss_stderr', 'f_offset_ss', 'f_offset_ss_stderr',
    'd_offset_ss', 'd_offset_ss_stderr',
    'Work_(pN*nm)', 'Work_(kB*T)',
    'delta Lc', 'total Lc', 'total W', 'total number of steps',
    'fit_status'
]


def make_failed_fit_row(status_msg):
    """Schema-conformant placeholder for a fit that raised an exception."""
    row = {key: np.nan for key in FIT_HEADER}
    row['fit_status'] = status_msg
    return row

def write_total_rows(filename_total_results, filename_i, steps_df, fit_df):
    """Concatenate step + fit results and append them to the total CSV,
    guaranteeing the output columns exactly match the file header."""
    expected = list(steps_df.columns) + FIT_HEADER

    if list(fit_df.columns) != FIT_HEADER:
        print(f"WARNING [{filename_i}]: fit columns mismatch, re-aligning: "
              f"{list(fit_df.columns)}")
        fit_df = fit_df.reindex(columns=FIT_HEADER)

    results_total_total = pd.concat([steps_df, fit_df], axis=1)

    if list(results_total_total.columns) != expected:
        print(f"WARNING [{filename_i}]: column mismatch after concat: "
              f"{list(results_total_total.columns)}")

    results_total_total.to_csv(filename_total_results, mode='a',
                               index=False, header=False)



"""define the functions of the subprocess processing the data"""


def show_h5_structure(file_path):
    file_h5 = lk.File(file_path)

    return file_h5


def read_in_data(file_num, Files, input_settings, input_format):
        # --- NEW: Check for Trap Position data ---
    use_trap_pos = False
    trap_pos_fd = None
    trap_pos_fd_ds = None
    if input_format['CSV'] == 1:
        df = pd.read_csv(Files[file_num])
        directory_i = Path(Files[file_num])
        filename_i = directory_i.name[:-4]
        # access the raw data
        Force = df.to_numpy()[:, 0]
        if input_format['length_measure'] == 1:
            Distance = df.to_numpy()[:, 1]
        else:
            Distance = df.to_numpy()[:, 1] / 1000
        # accessing the data frequency from user input
        Frequency_value = input_settings['data_frequency']  
        
        Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)
    
    else:
        with h5py.File(Files[file_num], "r") as f:
            directory_i = Path(Files[file_num])
            filename_i = directory_i.name[:-3]

            # access the raw data
            if input_format['HF'] == 1:
                if input_format['Trap'] == 1:
                    Force = f.get("Force HF/Force 1x")
                elif input_format['Trap'] == 0:
                    Force = f.get("Force HF/Force 2x")
                Distance = f.get("Distance/Piezo Distance")
                # accessing the data frequency from the h5 file
                Frequency_value = Force.attrs['Sample rate (Hz)']
                Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)


            elif input_format['Min step length'] == 1:
                if input_format['Trap'] == 1:
                    Force = f.get("Force HF/Force 1x")
                elif input_format['Trap'] == 0:
                    Force = f.get("Force HF/Force 2x")
                Distance = f.get("Distance/Piezo Distance")
                # accessing the data frequency from the h5 file
                Frequency_value = Force.attrs['Sample rate (Hz)']
                Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)

            elif input_format['LF'] == 1:
                if input_format['Trap'] == 1:
                    load_force = f.get("Force LF/Force 1x")
                    Force = load_force[:]['Value'][:]
                    try:
                        load_distance = f.get("Distance/Distance 1x")[:]
                    except:
                        load_distance = f.get("Distance/Distance 2")[:]
                    Distance = load_distance['Value'][:]
                elif input_format['Trap'] == 0:
                    load_force = f.get("Force LF/Force 2x")
                    Force = load_force[:]['Value'][:]
                    try:
                        load_distance = f.get("Distance/Distance 2x")[:]
                    except:
                        load_distance = f.get("Distance/Distance 1")[:]
                    Distance = load_distance['Value'][:]

                Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)

                # calculating the data frequency based on start- and end-time of the measurement
                size_F_LF = len(Force)
                stop_time_F_LF = load_force.attrs['Stop time (ns)']
                timestamp_F_LF = load_force.attrs['Start time (ns)']
                Frequency_value = size_F_LF / ((stop_time_F_LF - timestamp_F_LF) / 10**9)

            # VALIDATION: Ensure data arrays are not empty before returning
            if len(Force) == 0 or len(Distance) == 0:
                raise ValueError(f"No data found in {Files[file_num]}")


            if input_format['HF'] == 1 or input_format['Min step length'] == 1:
                try:
                    trap_path = "Trap position/1X"
                    if trap_path in f:
                        trap_raw = f.get(trap_path)
                        use_trap_pos = True

                        # Pass trap position through the SAME preprocessing pipeline as Distance
                        # This guarantees identical length, downsampling, and filtering
                        trap_pos_fd, _, trap_pos_fd_ds = preprocess_RAW(
                            Force, trap_raw, input_settings, input_format
                        )

                        print(f"Trap Position detected. Will use for curve splitting.")
                        print(f"  Force_Distance length: {len(Force_Distance)}, Trap FD length: {len(trap_pos_fd)}")
                    else:
                        print("No Trap Position data found. Using Distance for splitting.")
                except Exception as e:
                    print(f"Warning: Could not load Trap Position: {e}")

    return Force_Distance, Force_Distance_um, Frequency_value, filename_i, Force_Distance_ds, use_trap_pos, trap_pos_fd, trap_pos_fd_ds


def read_in_data_TOMATO(file_num, Files, input_settings, input_format):
    if input_format['CSV'] == 1:
        df = pd.read_csv(Files[file_num])
        directory_i = Path(Files[file_num])
        filename_i = directory_i.name[:-4]
        # access the raw data
        Force = df.to_numpy()[:, 0]
        if input_format['length_measure'] == 1:
            Distance = df.to_numpy()[:, 1]
        else:
            Distance = df.to_numpy()[:, 1] / 1000
        # accessing the data frequency from user input
        Frequency_value = input_settings['data_frequency']  
        
        Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)
    
    else:
        with h5py.File(Files[file_num], "r") as f:
            directory_i = Path(Files[file_num])
            filename_i = directory_i.name[:-3]

            # access the raw data
            if input_format['HF'] == 1:
                if input_format['Trap'] == 1:
                    Force = f.get("Force HF/Force 1x")
                elif input_format['Trap'] == 0:
                    Force = f.get("Force HF/Force 2x")
                Distance = f.get("Distance/Piezo Distance")
                # accessing the data frequency from the h5 file
                Frequency_value = Force.attrs['Sample rate (Hz)']
                Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)


            elif input_format['Min step length'] == 1:
                if input_format['Trap'] == 1:
                    Force = f.get("Force HF/Force 1x")
                elif input_format['Trap'] == 0:
                    Force = f.get("Force HF/Force 2x")
                Distance = f.get("Distance/Piezo Distance")
                # accessing the data frequency from the h5 file
                Frequency_value = Force.attrs['Sample rate (Hz)']
                Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)

            elif input_format['LF'] == 1:
                if input_format['Trap'] == 1:
                    load_force = f.get("Force LF/Force 1x")
                    Force = load_force[:]['Value'][:]
                    try:
                        load_distance = f.get("Distance/Distance 1x")[:]
                    except:
                        load_distance = f.get("Distance/Distance 2")[:]
                    Distance = load_distance['Value'][:]
                elif input_format['Trap'] == 0:
                    load_force = f.get("Force LF/Force 2x")
                    Force = load_force[:]['Value'][:]
                    try:
                        load_distance = f.get("Distance/Distance 2x")[:]
                    except:
                        load_distance = f.get("Distance/Distance 1")[:]
                    Distance = load_distance['Value'][:]

                Force_Distance, Force_Distance_um, Force_Distance_ds = preprocess_RAW(Force, Distance, input_settings, input_format)

                # calculating the data frequency based on start- and end-time of the measurement
                size_F_LF = len(Force)
                stop_time_F_LF = load_force.attrs['Stop time (ns)']
                timestamp_F_LF = load_force.attrs['Start time (ns)']
                Frequency_value = size_F_LF / ((stop_time_F_LF - timestamp_F_LF) / 10**9)

            

    return Force_Distance, Force_Distance_um, Frequency_value, filename_i, Force_Distance_ds

# open a folder containing raw data and lead through the analysis process
def start_subprocess(analysis_folder, timestamp, Files, input_settings, input_format, export_data, input_fitting, output_q):
    # create file to store total results
    if export_data['export_TOTAL'] == 1:
        filename_total_results = analysis_folder + '/total_results_' + timestamp + '.csv'

        with open(filename_total_results, 'w') as f:
            head = (
                'filename', 'orientation', 'Derivative of', 'step number',
                'F1', 'F2', 'Fc', 'step start', 'step end', 'step length'
            ) + tuple(FIT_HEADER)
            f.write(','.join(head))
            f.write('\n')

    # iterate through the files in the selected folder
    file_num = 0
    while file_num < len(Files):
        if file_num == 0:
            print('\nHard work ahead!\n')
            output_q.put('Hard work ahead!')

        # proceed differently with h5 and csv files
        Force_Distance, Force_Distance_um, Frequency_value, filename, Force_Distance_ds,  use_trap_pos, trap_pos_fd, trap_pos_fd_ds = read_in_data(file_num, Files, input_settings, input_format)


        # VALIDATION: Check if data was loaded successfully
        if len(Force_Distance) == 0:
            print(f"WARNING: File {filename} returned empty Force_Distance array! Skipping...")
            output_q.put(f'Error: File {filename} contains no valid data - skipped')
            file_num = file_num + 1
            continue


        num_curves = 1

        ###### Detect MultiFiles ######
        if input_format['MultiH5'] == 1:
            try:
                fw_curves, rv_curves, fw_curves_ds, rv_curves_ds = split_H5(
                                                                        Force_Distance, 
                                                                        Force_Distance_ds, 
                                                                        input_settings, 
                                                                        Frequency_value,
                                                                        use_trap_pos=use_trap_pos,           # New parameter
                                                                        trap_data=trap_pos_fd,               # New parameter (create this array)
                                                                        trap_ds_data=trap_pos_fd_ds          # New parameter (create this array)
                                                                    )
                                
                # Convert the returned object-arrays to Python lists of numpy arrays
                # Explicitly ensure each element is a numpy array
                fw_list = [np.asarray(item) for item in fw_curves] if hasattr(fw_curves, '__iter__') else [np.asarray(fw_curves)]
                rv_list = [np.asarray(item) for item in rv_curves] if hasattr(rv_curves, '__iter__') else [np.asarray(rv_curves)]
                fw_ds_list = [np.asarray(item) for item in fw_curves_ds] if hasattr(fw_curves_ds, '__iter__') else [np.asarray(fw_curves_ds)]
                rv_ds_list = [np.asarray(item) for item in rv_curves_ds] if hasattr(rv_curves_ds, '__iter__') else [np.asarray(rv_curves_ds)]

                num_fw = len(fw_list)
                num_rv = len(rv_list)

                if len(rv_list) == 0 or len(fw_list) == 0:
                    raise ValueError('No forward or no reverse curve found!')

                # Concatenate the lists
                curves = fw_list + rv_list
                curves_ds = fw_ds_list + rv_ds_list
                
                num_curves = len(curves)
                print(f"Success: Detected {num_curves} curves.")
                # Debug: check the type of the first element
                if num_curves > 0:
                    print(f"Type of curves[0]: {type(curves[0])}")
                    print(f"Shape of curves[0]: {curves[0].shape}")
                    print(f"Is curves[0] a list? {isinstance(curves[0], list)}")

            except Exception as e: 
                print(f"Error in MultiH5 processing: {e}")
                print('No Multi-File detected! Falling back to single curve.')
                curves = [Force_Distance]
                curves_ds = [Force_Distance_ds]
                num_curves = 1
        else:
            curves = [Force_Distance]
            curves_ds = [Force_Distance_ds]
            num_curves = 1
        print("number of curves is:")
        print(num_curves)
        print("type of curves")
        print(type(curves))
        

        for x in range(num_curves):
            # empty dataframe to store all step results of all curves in the folder
            total_results_steps = pd.DataFrame()

            # create dataframe to store all fitting parameters of all curves in the folder

            total_results_fit = pd.DataFrame(columns=FIT_HEADER)

            if num_curves == 1:
                filename_i = filename
            else:
                if x < num_fw:
                    suffix = 'fw_curve{num}'.format(num=x + 1)
                    filename_i = filename + '_' + suffix
                else:
                    suffix = 'rv_curve{num}'.format(num=x + 1 - num_fw)
                    filename_i = filename + '_' + suffix

            Force_Distance = curves[x][:, :2]
            print('################ FD', len(Force_Distance))
            Force_Distance_um = np.copy(Force_Distance)
            Force_Distance_um[:, 1] = Force_Distance_um[:, 1] / 1000
            Force_Distance_ds = curves_ds[x][:, :2]
        ###### Detect MultiFiles end ######

            orientation = "forward"
            # SAFE orientation check - prevent IndexError on empty arrays
            if len(Force_Distance) < 2:
                print(f"WARNING: Not enough data points ({len(Force_Distance)}) for orientation check in {filename_i}")
                output_q.put(f'Warning: Insufficient data points in {filename_i}')
            elif Force_Distance[0, 1] > Force_Distance[-1, 1]:  # reverse
                orientation = "reverse"
                Force_Distance = np.flipud(Force_Distance)
                Force_Distance_ds = np.flipud(Force_Distance_ds)
                Force_Distance_um = np.flipud(Force_Distance_um)

            # Export down sampled and smoothened FD values
            if export_data['export_SMOOTH'] == 1:
                save_to = analysis_folder + "/" + filename_i + "_smooth_" + timestamp + ".csv"
                with open(save_to, "w") as f:
                    np.savetxt(f, Force_Distance_um, delimiter=",")
                #np.savetxt(save_to, Force_Distance_um, delimiter=",")
            else:
                pass

            # trim data below specified force thresholds
            F_trimmed, PD_trimmed, F_low = trim_data(Force_Distance, input_settings['F_min'])
            common_steps = []
            print('#################### Trimmmed', len(F_trimmed))
            if not F_trimmed.size == 0 and not F_trimmed.size == 1:
                # create force and distance derivative of the pre-processed data to be able to identify steps
                derivative_array = create_derivative(input_settings, Frequency_value, F_trimmed, PD_trimmed, F_low)
                print('################### der array', len(derivative_array))

                # --- VALIDATION: skip curve if derivative array is invalid ---
                skip_curve = False
                if len(derivative_array) < 2:
                    print(f"WARNING: derivative_array too short ({len(derivative_array)}) for {filename_i}. Skipping curve.")
                    output_q.put(f'Warning: Skipped {filename_i} - insufficient data for derivative')
                    skip_curve = True
                elif len(derivative_array.shape) < 2:
                    print(f"WARNING: derivative_array is 1D for {filename_i}. Skipping curve.")
                    output_q.put(f'Warning: Skipped {filename_i} - invalid derivative array shape')
                    skip_curve = True

                if skip_curve:
                    # Write a dummy row to total_results so CSV export doesn't break
                    results_total_total = pd.concat([
                        pd.DataFrame({'filename': filename_i}, index=[0]),
                        pd.DataFrame(columns=FIT_HEADER)
                    ], axis=1)
                    if export_data['export_TOTAL'] == 1:
                        results_total_total.to_csv(filename_total_results, mode='a', index=False, header=False)
                    output_q.put(f'Done: Skipped {filename_i}')
                    continue
                # --- END VALIDATION ---

                """find steps based on force derivative"""
                filename_results = analysis_folder + "/" + filename_i + "_results_" + timestamp + ".csv"

                if input_format['Min_step_length'] == 1:
                    try:
                        results_F, PD_start_F = find_steps_F(
                            input_settings,
                            filename_i,
                            Force_Distance,
                            derivative_array,
                            orientation
                        )
                    except Exception as e:
                        print(f"Error in find_steps_F for {filename_i}: {e}")
                        traceback.print_exc()
                        results_F, PD_start_F = [], []
                else:
                    try:
                        results_F, PD_start_F = find_steps_F_old(
                            input_settings,
                            filename_i,
                            Force_Distance,
                            derivative_array,
                            orientation
                        )
                    except Exception as e:
                        print(f"Error in find_steps_F_old for {filename_i}: {e}")
                        traceback.print_exc()
                        results_F, PD_start_F = [], []

                results_F_list = list(results_F)

                if export_data['export_STEPS'] == 1:
                    steps_results_F = pd.DataFrame(results_F_list)
                    #with open(filename_results, 'a+') as f:
                        #f.write('\nSteps found by force derivative:\n')
                    steps_results_F.to_csv(filename_results, mode='a', index=False, header=True)
                else:
                    pass

                # except:
                #     results_F = []
                #     PD_start_F = []
                #     print("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Force steps')
                #     pass

                """find steps based on distance derivative"""

                if input_format['Min_step_length'] == 1:

                    try:
                        results_PD, PD_start_PD = find_steps_PD(
                            input_settings,
                            filename_i,
                            Force_Distance,
                            derivative_array,
                            orientation
                        )

                        results_PD_list = list(results_PD)

                        if export_data['export_STEPS'] == 1 and len(results_PD_list) > 0:
                            steps_results_PD = pd.DataFrame(results_PD_list)
                            steps_results_PD.to_csv(filename_results, mode='a', index=False, header=True)

                    except Exception as e:
                        results_PD = []
                        PD_start_PD = []
                        err_PD = str("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Distance steps')
                        print(err_PD)
                        print(f"Detailed error: {e}")
                        traceback.print_exc()

                else:
                    try:
                        results_PD, PD_start_PD = find_steps_PD_old(
                            input_settings,
                            filename_i,
                            Force_Distance,
                            derivative_array,
                            orientation
                        )

                        results_PD_list = list(results_PD)

                        if export_data['export_STEPS'] == 1 and len(results_PD_list) > 0:
                            steps_results_PD = pd.DataFrame(results_PD_list)
                            steps_results_PD.to_csv(filename_results, mode='a', index=False, header=True)

                    except Exception as e:
                        results_PD = []
                        PD_start_PD = []
                        err_PD = str("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Distance steps')
                        print(err_PD)
                        print(f"Detailed error: {e}")
                        traceback.print_exc()

                # save plot with FD-curve, derivatives and found steps
                # save plot with FD-curve, derivatives and found steps
                try:
                    save_figure(
                        export_data['export_PLOT'],
                        export_data, 
                        timestamp,
                        filename_i,
                        analysis_folder,
                        Force_Distance,
                        derivative_array,
                        F_trimmed,
                        PD_trimmed,
                        PD_start_F if isinstance(PD_start_F, list) else [],
                        PD_start_PD if isinstance(PD_start_PD, list) else []
                    )
                except Exception as e:
                    print(f"Warning: Could not save figure for {filename_i}: {e}")
                    traceback.print_exc()
                    pass



                # when steps are found by force AND distance derivative, they are considered common steps
                #common_steps = []
                common_steps = []
                try:
                    if len(results_F_list) > 0 and len(results_PD_list) > 0:
                        common_steps = find_common_steps(results_F_list, results_PD_list)
                    else:
                        print(f"No steps found for {filename_i} (F: {len(results_F_list)}, PD: {len(results_PD_list)})")

                    print("common steps are:")
                    print(common_steps)

                    common_steps_results = [{'filename': filename_i, 'orientation': orientation, 'Derivative of': '', 'step #': 0, 'F1': '', 'F2': '', 'Fc': '', 'step start': '', 'step end': '', 'step length': ''}]
                except Exception as e:
                    err_FCS = str("Error in finding common steps for " + str(filename_i) + ': ' + str(e))
                    output_q.put(err_FCS)
                    traceback.print_exc()

                # append common steps to the 'step 0'
                # append common steps to the 'step 0'
                if common_steps:
                    for cs in range(len(common_steps)):
                        common_steps_results.append(common_steps[cs])

                    # convert common steps to dataframe for export
                    common_steps_results = pd.DataFrame(common_steps_results)

                    # export the steps into the results for ONLY this file
                    #with open(filename_results, 'a+') as f:
                        #f.write('\nCommon steps:\n')
                    common_steps_results.to_csv(filename_results, mode='a', index=False, header=True)

                    # put common steps into a total_results dataframe so all steps from all files of the analysed folder can be exported together
                    total_results_steps = pd.concat([total_results_steps, common_steps_results], ignore_index=True, sort=False)

                else:
                    common_steps_results = pd.DataFrame({'filename': filename_i, 'orientation': orientation, 'Derivative of': '', 'step #': 0, 'F1': '', 'F2': '', 'Fc': '', 'step start': '', 'step end': '', 'step length': ''}, index=[0])
                    total_results_steps = pd.concat([total_results_steps, common_steps_results], ignore_index=True, sort=False)

                '''if common steps were found, try to fit FD-Curve'''

                if export_data['export_FIT'] == 1:
                    try:
                        export_fit = []
                        fit = []
                        start_force_ss = []
                        start_distance_ss = []
                        integral_ss_fit_start = []
                        integral_ss_fit_end = []

                        # try to fit all parts of curve based on the common steps
                        try:
                                    ###### Reverse fitting ######
                            if input_format['reverse_fitting'] == 1:
                                try:
                                    export_fit_ds_FU, area_ds, step_start, model_FU_obj = fitting_FU(
                                        filename_i,
                                        input_settings,
                                        export_data,
                                        input_fitting,
                                        float(common_steps[0]['step start']),
                                        float(common_steps[-1]['step end']),
                                        max(derivative_array[:, 1]),
                                        Force_Distance_ds,
                                        derivative_array,
                                        F_low,
                                        0
                                    )

                                    fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ds, area_ds = fitting_FU_ss(
                                            filename_i,
                                            input_settings,
                                            export_data,
                                            input_fitting,
                                            float(0),
                                            float(common_steps[0]['step start']),
                                            Force_Distance_ds,
                                            1,
                                            1,
                                            derivative_array,
                                            F_low,
                                            0
                                        )
                                    
                                    
                                    fit.append(fit_ss)
                                    start_force_ss.append(f_fitting_region_ss)
                                    start_distance_ss.append(d_fitting_region_ss)
                                    
                                   


                                except Exception as e:
                                    export_fit.append(make_failed_fit_row(f'error: {type(e).__name__}: {e}'))
                                    
                                    print(f"Error: {e}")
                                    traceback.print_exc()
                                    print('Something went wrong with reverse fitting!')
                                    
                            else:

                                # fit part between start of the FD-cure up to the first common step
                                export_fit_ds, area_ds, step_start, main_ds_model = fitting_ds(
                                    filename_i,
                                    input_settings,
                                    export_data,
                                    input_fitting,
                                    float(common_steps[0]['step start']),
                                    Force_Distance_ds,
                                    derivative_array,
                                    F_low,
                                    0
                                )





                            export_fit.append(export_fit_ds)

                            # fit parts after steps, when more than one common step was found, there are multiple parts to fit
                            if len(common_steps) > 1:
                                #print(common_steps)
                                print("length of common steps is "+str(len(common_steps)))
                                for n in range(0, len(common_steps) - 1):
                                    # try to fit each part of the curve, if one of the parts can not be fitted, still try to fit the others
                                    try:
                                        
                                        fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss, area_ss_fit_start, area_ss_fit_end = fitting_ss(
                                            filename_i,
                                            input_settings,
                                            export_data,
                                            input_fitting,
                                            float(common_steps[n]['step end']),
                                            float(common_steps[n + 1]['step start']),
                                            Force_Distance_ds,
                                            1,
                                            1,
                                            derivative_array,
                                            F_low,
                                            0
                                        )
                                        #print(str(n))
                                        fit.append(fit_ss)
                                        start_force_ss.append(f_fitting_region_ss)
                                        start_distance_ss.append(d_fitting_region_ss)
                                        export_fit.append(export_fit_ss)
                                        integral_ss_fit_start.append(area_ss_fit_start)
                                        integral_ss_fit_end.append(area_ss_fit_end)

                                    except Exception as e:
                                        export_fit.append(make_failed_fit_row(f'error: {type(e).__name__}: {e}'))
                                        print("something went wrong with the middle part of ss fitting")
                                        print(f"Error: {e}")
                                        traceback.print_exc()
                                        pass

                            # fit the last part of the curve
                            try:
                                fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss, area_ss_fit_start, area_ss_fit_end = fitting_ss(
                                    filename_i,
                                    input_settings,
                                    export_data,
                                    input_fitting,
                                    float(common_steps[len(common_steps) - 1]['step end']),
                                    max(derivative_array[:, 1]),
                                    Force_Distance_ds,
                                    1,
                                    1,
                                    derivative_array,
                                    F_low,
                                    0
                                )

                                fit.append(fit_ss)
                                start_force_ss.append(f_fitting_region_ss)
                                start_distance_ss.append(d_fitting_region_ss)
                                export_fit.append(export_fit_ss)
                                integral_ss_fit_start.append(area_ss_fit_start)
                                integral_ss_fit_end.append(area_ss_fit_end)

                            except Exception as e:
                                export_fit.append(make_failed_fit_row(f'error: {type(e).__name__}: {e}'))
                                print("something went wrong with the last part of ss fitting")
                                print(f"Error: {e}")
                                traceback.print_exc()
                                pass

                            '''from the fits, work put into the system is calculated'''
                            if common_steps:
                                work_per_step = [0]  # in pN*nm
                                kT_per_step = [0]    # in kT

                                work_first_step, kT_1 = calc_integral(
                                    area_ds,
                                    integral_ss_fit_start[0],
                                    common_steps[0]['step start'],
                                    common_steps[0]['step end'],
                                    common_steps[0]['F1'],
                                    common_steps[0]['F2']
                                )

                                #print("Work of first step: " + str(work_first_step))
                                work_per_step.append(work_first_step)
                                kT_per_step.append(kT_1)

                                if len(common_steps) > 1:
                                    for n in range(0, len(common_steps) - 1):
                                        work_step_n, kT_n = calc_integral(
                                            integral_ss_fit_end[n],
                                            integral_ss_fit_start[n + 1],
                                            common_steps[n + 1]['step start'],
                                            common_steps[n + 1]['step end'],
                                            common_steps[n + 1]['F1'],
                                            common_steps[n + 1]['F2']
                                        )

                                        work_per_step.append(work_step_n)
                                        kT_per_step.append(kT_n)
                                
                                j = 0
                                for dict in export_fit:
                                    dict["Work_(pN*nm)"] = work_per_step[j]
                                    dict["Work_(kB*T)"] = kT_per_step[j]
                                    j += 1
                                
                        # if no step was found, the common step index 0 is not available and will raise an IndexError.
                        # So in this case the fit will be performed for the whole curve from beginning to end.
                        except IndexError:
                            if not common_steps:
                                export_fit_ds, area_ds, step_start, main_ds_model = fitting_ds(
                                    filename_i,
                                    input_settings,
                                    export_data,
                                    input_fitting,
                                    derivative_array[-1, 1],
                                    Force_Distance_ds,
                                    derivative_array,
                                    F_low,
                                    0
                                )

                                export_fit.append(export_fit_ds)
                                print("no common steps found")
                        
                        # Build the fit DataFrame on the canonical schema:
                        # unknown keys are dropped, missing keys become NaN.
                        export_fit_df = pd.DataFrame(export_fit or [],
                                                     columns=FIT_HEADER)
                        export_fit_df = export_fit_df.reindex(columns=FIT_HEADER)

                        # Use pd.concat to append the data
                        total_results_fit = pd.concat([total_results_fit, export_fit_df], ignore_index=True, sort=False)


                        #total_results_fit = total_results_fit.append(export_fit, ignore_index=True, sort=False)
                       
                        # create a plot for the fitted curve
                        try:
                            plot_fit(fit, 
                                    start_force_ss, 
                                    start_distance_ss,
                                    Force_Distance, 
                                    Force_Distance_ds, 
                                    analysis_folder, 
                                    filename_i, 
                                    timestamp,
                                    export_data=export_data, 
                                    model_FU=model_FU_obj if 'model_FU_obj' in locals() else None,
                                    model_ds_final=main_ds_model if 'main_ds_model' in locals() else None
                                    )
                        except: 
                            pass
                    except Exception as e:
                        print(f"Error: {e}")
                        traceback.print_exc()
                        print('Something went wrong with fitting')
                        
                        pass
                
                #print("total results steps are:")
                #print(total_results_steps)



                #print("total results fits are:")
                #print(total_results_fit)

                # Remove the last three columns
                #total_results_fit = total_results_fit.iloc[:, :-3]

            
            
                """start of tab shift"""
            if not common_steps:
                print("no common steps found")
                common_steps = []
            if len(common_steps)>0 and len(F_trimmed) >0:
                #print(common_steps)

                if len(total_results_fit) > 0:
                    # Calculate delta LC
                    total_results_fit['delta Lc'] = total_results_fit['Lc_ss'].diff().fillna("#N/A")

                    # Initialize total Lc and total W columns with #N/A
                    total_results_fit['total Lc'] = "#N/A"
                    total_results_fit['total W'] = "#N/A"
                    total_results_fit['total number of steps'] = "#N/A"

                    # Set total Lc for the last row
                    try:
                        total_results_fit.loc[total_results_fit.index[-1], 'total Lc'] = total_results_fit['Lc_ss'].iloc[-1]
                    except (IndexError, KeyError):
                        pass

                    # Set total W for the last row
                    try:
                        total_results_fit.loc[total_results_fit.index[-1], 'total W'] = total_results_fit['Work_(kB*T)'].sum()
                    except (IndexError, KeyError):
                        pass

                    # Set total Lc for the last row
                    try:
                        total_results_fit.loc[total_results_fit.index[-1], 'total number of steps'] = total_results_steps['step #'].iloc[-1]
                    except (IndexError, KeyError):
                        pass


                write_total_rows(filename_total_results, filename_i,
                                 total_results_steps, total_results_fit)

                print('done', x + 1, 'curves from', len(curves))
                out_progress = str('File ' + str(file_num + 1) + ': Done ' + str(x + 1) + ' curves from ' + str(len(curves)))
                output_q.put(out_progress)

                print(filename_i)
                output_q.put(filename_i)
            else:
                write_total_rows(filename_total_results, filename_i,
                                                total_results_steps, total_results_fit)

                print('This curve was below the Force threshold and could not be processed!\nPlease check if the correct trap was selected.')
                output_q.put('This curve was below the Force threshold and could not be processed!\nPlease check if the correct trap was selected.')
            
                """end of tab shift"""




        if file_num == int(len(Files) / 2):
            print('\nHalf way there!\n')
            output_q.put('Half way there!')
            print()
        elif file_num == len(Files) - 1:
            print('\nAlmost there!\n')
            output_q.put('Almost there!')

        file_num = file_num + 1
        print('done', file_num, 'from', len(Files))
        out_progress = str('Done ' + str(file_num) + ' files from ' + str(len(Files)))
        output_q.put(out_progress)

            # --- MEMORY CLEANUP START ---
        try:
            del Force_Distance, Force_Distance_um, Force_Distance_ds, \
                derivative_array, F_trimmed, PD_trimmed, common_steps, \
                results_F_list, results_PD_list, export_fit, fit, \
                start_force_ss, start_distance_ss
        except NameError:
            # Some names were never created this iteration
            # (skipped curve / below force threshold)
            pass

        gc.collect()
        # --- MEMORY CLEANUP END ---

        print(filename_i)
        output_q.put(filename_i)

    print('Analysis finished! \nProgram can be closed.')
    output_q.put('Analysis finished! \nProgram can be closed.')
