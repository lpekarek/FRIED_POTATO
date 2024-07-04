"""FRIED POTATO -- 2024 -- Version 1"""

import pandas as pd
import h5py
import numpy as np
from pathlib import Path
import traceback
import lumicks.pylake as lk
from scipy import signal

# relative imports
from FRIED_POTATO_fitting import fitting_ds, fitting_ss, plot_fit, fitting_FU, fitting_FU_ss
from FRIED_POTATO_preprocessing import trim_data, create_derivative
from FRIED_POTATO_find_steps import find_steps_F, find_steps_PD, find_common_steps, calc_integral, save_figure
from FRIED_POTATO_processMultiH5 import split_H5

"""define the functions of the subprocess processing the data"""


class force_ramp_data():
    def __init__(self, file_num, Files, input_text, input_checkbox):
        self.directory_i = Path(Files[file_num])
        self.filename_i = self.directory_i.stem
        self.force_distance = np.empty((0, 2)) # initialize empty array with 2 columns
        self.output = []

        if input_checkbox['CSV'] == 1:
            self.df = pd.read_csv(Files[file_num])

            # access the raw data
            force = self.df.to_numpy()[:, 0]
            if input_checkbox['nm'] == 1:
                distance = self.df.to_numpy()[:, 1]
            else:
                distance = self.df.to_numpy()[:, 1] / 1000
            # accessing the data frequency from user input
            frequency_value = input_text['frequency']  

        else:
            h5_file = lk.File(Files[file_num])

            try:
                if input_checkbox['1x']:
                    force = h5_file['Force HF']['Force 1x'].data
                elif input_checkbox['2x']:
                    force = h5_file['Force HF']['Force 2x'].data
                distance = h5_file['Distance']['Piezo Distance'].data
            
            except:
                if input_checkbox['1x']:
                    force = h5_file['Force LF']['Force 1x'].data
                    distance = h5_file['Distance']['Distance 1x'].data
                elif input_checkbox['2x']:
                    force = h5_file['Force LF']['Force 2x'].data
                    distance = h5_file['Distance']['Distance 2x'].data

            #     # calculating the data frequency based on start- and end-time of the measurement
            #     size_F_LF = len(force)
            #     stop_time_F_LF = load_force.attrs['Stop time (ns)']
            #     timestamp_F_LF = load_force.attrs['Start time (ns)']
            #     Frequency_value = size_F_LF / ((stop_time_F_LF - timestamp_F_LF) / 10**9)
            
            # # accessing the data frequency from the h5 file
            # frequency_value = force.attrs['Sample rate (Hz)']
        self.force_distance = np.column_stack((force, distance))
        self.output.append(f'Data points before preprocessing: {len(self.force_distance):,}')
        self.preprocess_RAW(input_text, input_checkbox)
        self.output.append(f'Data points after preprocessing: {len(self.force_distance):,}')


    def preprocess_RAW(self, input_text, input_checkbox):
        if input_checkbox['aug'] == 1 and input_text['augment_factor'] != '':
            # Augment data
            new_f = []
            new_d = []
            factor = int(input_text['augment_factor'])

            for line in range(len(self.force_distance) - 1):
                f = self.force_distance[line, 0]
                d = self.force_distance[line, 1]
                new_f.append(f)
                new_d.append(d)
                f_next = self.force_distance[line + 1, 0]
                d_next = self.force_distance[line + 1, 1]
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

            self.force_distance[:, 0] = np.array(new_f)
            self.force_distance[:, 1] = np.array(new_d)


        if input_checkbox['prepro'] == 1:
            # Downsample
            self.force_distance = self.force_distance[::int(input_text['downsample'])]

            # Filter
            b, a = signal.butter(int(input_text['filter_degree']), float(input_text['filter_cut_off']))
            self.force_distance = signal.filtfilt(b, a, self.force_distance, axis=0)

# def read_in_data(file_num, Files, input_text, input_checkbox):
#     if input_checkbox['CSV'] == 1:
#         df = pd.read_csv(Files[file_num])
#         directory_i = Path(Files[file_num])
#         filename_i = directory_i.name[:-4]
#         # access the raw data
#         force = df.to_numpy()[:, 0]
#         if input_checkbox['length_measure'] == 1:
#             distance = df.to_numpy()[:, 1]
#         else:
#             distance = df.to_numpy()[:, 1] / 1000
#         # accessing the data frequency from user input
#         Frequency_value = input_text['frequency']  
#         force_distance, force_distance_um = preprocess_RAW(force, distance, input_text, input_checkbox)

#     else:
#         h5_file = lk.File(Files[file_num]
# )
#         with h5py.File(Files[file_num], "r") as f:
#             directory_i = Path(Files[file_num])
#             filename_i = directory_i.name[:-3]

#             # access the raw data
#             if input_checkbox['HF'] == 1:
#                 if input_checkbox['1x'] == 1:
#                     force = f.get("force HF/force 1x")
#                 elif input_checkbox['2x'] == 1:
#                     force = f.get("force HF/force 2x")
#                 distance = f.get("distance/Piezo distance")
#                 # accessing the data frequency from the h5 file
#                 Frequency_value = force.attrs['Sample rate (Hz)']
#                 force_distance, force_distance_um = preprocess_RAW(force, distance, input_text, input_checkbox)

#             elif input_checkbox['LF'] == 1:
#                 if input_checkbox['1x'] == 1:
#                     load_force = f.get("force LF/force 1x")
#                     force = load_force[:]['Value'][:]
#                     try:
#                         load_distance = f.get("distance/distance 1x")[:]
#                     except:
#                         load_distance = f.get("distance/distance 2")[:]
#                     distance = load_distance['Value'][:]
#                 elif input_checkbox['2x'] == 1:
#                     load_force = f.get("force LF/force 2x")
#                     force = load_force[:]['Value'][:]
#                     try:
#                         load_distance = f.get("distance/distance 2x")[:]
#                     except:
#                         load_distance = f.get("distance/distance 1")[:]
#                     distance = load_distance['Value'][:]

#                 force_distance, force_distance_um = preprocess_RAW(force, distance, input_text, input_checkbox)

#                 # calculating the data frequency based on start- and end-time of the measurement
#                 size_F_LF = len(force)
#                 stop_time_F_LF = load_force.attrs['Stop time (ns)']
#                 timestamp_F_LF = load_force.attrs['Start time (ns)']
#                 Frequency_value = size_F_LF / ((stop_time_F_LF - timestamp_F_LF) / 10**9)

#     return force_distance, force_distance_um, Frequency_value, filename_i


# open a folder containing raw data and lead through the analysis process
def start_subprocess(analysis_folder, timestamp, Files, input_settings, input_format, export_data, input_fitting, output_q):
    # create file to store total results
    if export_data['export_TOTAL'] == 1:
        filename_total_results = analysis_folder + '/total_results_' + timestamp + '.csv'

        with open(filename_total_results, 'w') as f:
            #f.write('>Common steps from all curves of the folder:\n') #messes up with reading this file later
            head = (
                'filename',
                'orientation',
                'Derivative of',
                'step number',
                'F1',
                'F2',
                'Fc',
                'step start',
                'step end',
                'step length',
                'filename',
                'model',
                'model_ss',
                'log_likelihood',
                'Lc_ds',
                'Lp_ds',
                'Lp_ds_stderr',
                'St_ds',
                'Lc_ss',
                'Lc_ss_stderr',
                'Lp_ss',
                'St_ss',
                'f_offset',
                'd_offset',
                'Work_(pN*nm)',
                'Work_(kB*T)', 
                "delta Lc",
                "total Lc",
                "total W"
            )
            f.write(','.join(head))
            f.write('\n')

    # iterate through the files in the selected folder
    file_num = 0
    while file_num < len(Files):
        if file_num == 0:
            print('\nHard work ahead!\n')
            output_q.put('Hard work ahead!')

        # proceed differently with h5 and csv files
        force_distance, force_distance_um, Frequency_value, filename = read_in_data(file_num, Files, input_settings, input_format)

        num_curves = 1

        ###### Detect MultiFiles ######
        if input_format['MultiH5'] == 1:
            try:
                fw_curves, rv_curves = split_H5(force_distance, input_settings, Frequency_value)
                num_fw = len(fw_curves)
                fw_curves.extend(rv_curves)
                curves = fw_curves
            except:
                print('No Multi-File detected!')
                curves = [force_distance]
        else:
            curves = [force_distance]

        num_curves = len(curves)

        for x in range(num_curves):
            # empty dataframe to store all step results of all curves in the folder
            total_results_steps = pd.DataFrame()

            # create dataframe to store all fitting parameters of all curves in the folder
            header_fit = [
                "filename",
                "model",
                'model_ss',
                "log_likelihood",
                "Lc_ds",
                "Lp_ds",
                "Lp_ds_stderr",
                "St_ds",
                "Lc_ss",
                "Lc_ss_stderr",
                "Lp_ss",
                "St_ss",
                "f_offset",
                "d_offset",
                "Work_(pN*nm)",
                "Work_(kB*T)",
                "delta Lc",
                "total Lc",
                "total W"
            ]

            total_results_fit = pd.DataFrame(columns=header_fit)

            if num_curves == 1:
                filename_i = filename
            else:
                if x < num_fw:
                    suffix = 'fw_curve{num}'.format(num=x + 1)
                    filename_i = filename + '_' + suffix
                else:
                    suffix = 'rv_curve{num}'.format(num=x + 1 - num_fw)
                    filename_i = filename + '_' + suffix

            force_distance = curves[x][:, :2]
            print('################ FD', len(force_distance))
            force_distance_um = np.copy(force_distance)
            force_distance_um[:, 1] = force_distance_um[:, 1] / 1000
        ###### Detect MultiFiles end ######

            orientation = "forward"
            if force_distance[0, 1] > force_distance[-1, 1]:  # reverse
                orientation = "reverse"
                force_distance = np.flipud(force_distance)
                force_distance_um = np.flipud(force_distance_um)

            # Export down sampled and smoothened FD values
            if export_data['export_SMOOTH'] == 1:
                save_to = analysis_folder + "/" + filename_i + "_smooth_" + timestamp + ".csv"
                np.savetxt(save_to, force_distance_um, delimiter=",")
            else:
                pass

            # trim data below specified force thresholds
            F_trimmed, PD_trimmed, F_low = trim_data(force_distance, input_settings['F_min'])
            print('#################### Trimmmed', len(F_trimmed))
            if not F_trimmed.size == 0:
                # create force and distance derivative of the pre-processed data to be able to identify steps
                derivative_array = create_derivative(input_settings, Frequency_value, F_trimmed, PD_trimmed, F_low)
                print('################### der array', len(derivative_array))
                """find steps based on force derivative"""
                filename_results = analysis_folder + "/" + filename_i + "_results_" + timestamp + ".csv"

                # try:
                results_F, PD_start_F = find_steps_F(
                    input_settings,
                    filename_i,
                    force_distance,
                    derivative_array,
                    orientation
                )

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
                #     print("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding force steps')
                #     pass

                """find steps based on distance derivative"""

                try:
                    results_PD, PD_start_PD = find_steps_PD(
                        input_settings,
                        filename_i,
                        force_distance,
                        derivative_array,
                        orientation
                    )

                    results_PD_list = list(results_PD)

                    if export_data['export_STEPS'] == 1:
                        steps_results_PD = pd.DataFrame(results_PD_list)
                        #with open(filename_results, 'a+') as f:
                            #f.write('\nSteps found by distance derivative:\n')
                        steps_results_PD.to_csv(filename_results, mode='a', index=False, header=True)

                except:
                    results_PD = []
                    PD_start_PD = []
                    err_PD = str("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding distance steps')
                    print(err_PD)
                    pass

                # save plot with FD-curve, derivatives and found steps
                save_figure(
                    export_data['export_PLOT'],
                    timestamp,
                    filename_i,
                    analysis_folder,
                    force_distance,
                    derivative_array,
                    F_trimmed,
                    PD_trimmed,
                    PD_start_F,
                    PD_start_PD
                )

                # when steps are found by force AND distance derivative, they are considered common steps
                common_steps = []
                try:
                    common_steps = find_common_steps(results_F_list, results_PD_list)
                    # to match with the fitting rows (always one more than steps) put a 'step 0' as first line
                    common_steps_results = [{'filename': filename_i, 'orientation': orientation, 'Derivative of': '', 'step #': 0, 'F1': '', 'F2': '', 'Fc': '', 'step start': '', 'step end': '', 'step length': ''}]
                except:
                    err_FCS = str("Error in finding common steps" + str(filename_i) + '\n' 'There was an error in finding common steps')
                    output_q.put(err_FCS)
                    pass

                # append common steps to the 'step 0'
                if common_steps:
                    for x in range(len(common_steps)):
                        common_steps_results.append(common_steps[x])

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
                empty = {
                    'filename': filename_i,
                    'model': 'None',
                    'log_likelihood': 'None',
                    'Lc_ds': 'None',
                    'Lp_ds': 'None',
                    'Lp_ds_stderr': 'None',
                    'St_ds': 'None',
                    'f_offset': 'None',
                    'd_offset': 'None'
                }

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
                                    export_fit_ds_FU, area_ds, step_start = fitting_FU(
                                        filename_i,
                                        input_settings,
                                        export_data,
                                        input_fitting,
                                        float(common_steps[0]['step start']),
                                        float(common_steps[-1]['step end']),
                                        max(derivative_array[:, 1]),
                                        force_distance,
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
                                            force_distance,
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
                                    export_fit.append(empty)
                                    
                                    print(f"Error: {e}")
                                    traceback.print_exc()
                                    print('Something went wrong with reverse fitting!')
                                    
                            else:

                                # fit part between start of the FD-cure up to the first common step
                                export_fit_ds, area_ds, step_start = fitting_ds(
                                    filename_i,
                                    input_settings,
                                    export_data,
                                    input_fitting,
                                    float(common_steps[0]['step start']),
                                    force_distance,
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
                                            force_distance,
                                            1,
                                            1,
                                            derivative_array,
                                            F_low,
                                            0
                                        )
                                        print(str(n))
                                        fit.append(fit_ss)
                                        start_force_ss.append(f_fitting_region_ss)
                                        start_distance_ss.append(d_fitting_region_ss)
                                        export_fit.append(export_fit_ss)
                                        integral_ss_fit_start.append(area_ss_fit_start)
                                        integral_ss_fit_end.append(area_ss_fit_end)

                                    except Exception as e:
                                        export_fit.append(empty)
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
                                    force_distance,
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
                                export_fit.append(empty)
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

                                print("Work of first step: " + str(work_first_step))
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
                                export_fit_ds, area_ds, step_start = fitting_ds(
                                    filename_i,
                                    input_settings,
                                    export_data,
                                    input_fitting,
                                    derivative_array[-1, 1],
                                    force_distance,
                                    derivative_array,
                                    F_low,
                                    0
                                )

                                export_fit.append(export_fit_ds)
                                print("no common steps found")
                        
                        # If export_fit is a list of dictionaries, convert it to a DataFrame
                        if isinstance(export_fit, list):
                            export_fit_df = pd.DataFrame(export_fit)
                        else:
                            export_fit_df = export_fit

                        # Use pd.concat to append the data
                        total_results_fit = pd.concat([total_results_fit, export_fit_df], ignore_index=True, sort=False)


                        #total_results_fit = total_results_fit.append(export_fit, ignore_index=True, sort=False)
                       
                        # create a plot for the fitted curve
                        plot_fit(fit, start_force_ss, start_distance_ss, force_distance, analysis_folder, filename_i, timestamp)
                        
                    except Exception as e:
                        print(f"Error: {e}")
                        traceback.print_exc()
                        print('Something went wrong with fitting')
                        
                        pass
                
                print("total results steps are:")
                print(total_results_steps)



                print("total results fits are:")
                print(total_results_fit)

                # Remove the last three columns
                #total_results_fit = total_results_fit.iloc[:, :-3]

                # Calculate delta LC
                total_results_fit['delta Lc'] = total_results_fit['Lc_ss'].diff().fillna("#N/A")

                # Initialize total Lc and total W columns with #N/A
                total_results_fit['total Lc'] = "#N/A"
                total_results_fit['total W'] = "#N/A"

                # Set total Lc for the last row
                total_results_fit.loc[total_results_fit.index[-1], 'total Lc'] = total_results_fit['Lc_ss'].iloc[-1]

                # Set total W for the last row
                total_results_fit.loc[total_results_fit.index[-1], 'total W'] = total_results_fit['Work_(kB*T)'].sum()

                # Find the index of the "Work_(kB*T)" column
                insert_pos = total_results_fit.columns.get_loc('Work_(kB*T)') + 1

                # Insert the new columns after "Work_(kB*T)"
                #total_results_fit = pd.concat([total_results_fit.iloc[:, :insert_pos], total_results_fit[['delta Lc', 'total Lc', 'total W']], total_results_fit.iloc[:, insert_pos:]], axis=1)


                print("total results fits are:")
                print(total_results_fit)



                results_total_total = pd.concat([total_results_steps, total_results_fit], axis=1)
                results_total_total.to_csv((filename_total_results), mode='a', index=False, header=False)

                print('done', x + 1, 'curves from', len(curves))
                out_progress = str('File ' + str(file_num + 1) + ': Done ' + str(x + 1) + ' curves from ' + str(len(curves)))
                output_q.put(out_progress)

                print(filename_i)
                output_q.put(filename_i)
            else:
                print('This curve was below the force threshold and could not be processed!\nPlease check if the correct trap was selected.')
                output_q.put('This curve was below the force threshold and could not be processed!\nPlease check if the correct trap was selected.')


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

        print(filename_i)
        output_q.put(filename_i)

    print('Analysis finished! \nProgram can be closed.')
    output_q.put('Analysis finished! \nProgram can be closed.')
