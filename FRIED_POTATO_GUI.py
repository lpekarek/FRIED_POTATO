
""" 
    FRIED POTATO -- 2024 -- Version 1
    Developed and maintained by Lukáš Pekárek and Stefan Buck

    Force-Ramp Improved EDition of Practical Optical Tweezers Analysis TOol

    Please also refer to the original POTATO publication: https://doi.org/10.1101/2021.11.11.468103

    FRIED POTATO was redesigned from the original POTATO to handle a broader range of use-cases.

    processes Force-Distance Optical Tweezers data in an automated way, to find unfolding events
    The tool is developed to handle h5 raw data, produced from the C-Trap OT machine from Lumicks,
    as well as any other FD data prepared in a csv file (2 columns: Force(pN) - Distance(um))
    Furthermore the script can analyse single constant force files
    The parameters can be changed in the GUI before each run.
    Alternatively they can be changed permanently in the POTATO_config file
"""

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from tkinter import ttk
import lumicks.pylake as lk
from PIL import ImageTk, Image
import pandas as pd
import numpy as np
import os
import glob
import time
import multiprocessing as mp
import json

# relative imports
from FRIED_POTATO_ForceRamp import start_subprocess, read_in_data
from FRIED_POTATO_preprocessing import create_derivative
from FRIED_POTATO_config import default_values_HF, default_values_LF, default_values_CSV, default_values_FIT, default_values_constantF
from FRIED_POTATO_constantF import get_constantF, display_constantF, fit_constantF
from FRIED_POTATO_fitting import fitting_ds, fitting_ss
from FRIED_POTATO_find_steps import calc_integral

# To avoid blurry GUI - DPI scaling
import ctypes
awareness = ctypes.c_int()
errorCode = ctypes.windll.shcore.GetProcessDpiAwareness(0, ctypes.byref(awareness))
errorCode = ctypes.windll.shcore.SetProcessDpiAwareness(1)


"""define the functions used in the GUI"""


# get settings, get folder directory, create analysis results folder
# def start_analysis():
#     global p0
#     global analysis_folder

#     # check user input
#     input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()

#     # ask wich directory should be analysed
#     folder = tk.filedialog.askdirectory()
#     root.title('FRIED POTATO -- ' + str(folder))

#     # decide which input format was choosen
#     if input_format['CSV'] == 1:
#         folder_path = str(folder + "/*.csv")
#     else:
#         folder_path = str(folder + "/*.h5")

#     Files = glob.glob(folder_path)

#     # print number of files to analyse, if no files found give an error
#     print('Files to analyse', len(Files))
#     output_window.insert("end", 'Files to analyse: ' + str(len(Files)) + "\n")
#     output_window.see("end")
#     if not len(Files) == 0:
#         output_window.insert("end", 'Analysis in progress. Please do not close the program! \n')

#         # print starting time of the analysis
#         timestamp = time.strftime("%Y%m%d-%H%M%S")
#         print("Timestamp: " + timestamp)
#         output_window.insert("end", 'Start of analysis: ' + str(timestamp) + "\n")
#         output_window.see("end")

#         # create a folder for the analysis results
#         analysis_folder = str(folder + '/Analysis_' + timestamp)
#         os.mkdir(analysis_folder)

#         # export configuration file with used parameters
#         export_settings(analysis_folder, timestamp, input_settings, input_fitting)

#         # start analysis in a new process
#         p0 = mp.Process(target=start_subprocess, name='Process-0', args=(
#             analysis_folder,
#             timestamp,
#             Files,
#             input_settings,
#             input_format,
#             export_data,
#             input_fitting,
#             output_q,
#         ))

#         p0.daemon = True
#         p0.start()

#     else:
#         output_window.insert("end", 'No file of the selected data type in the folder! \n')
#         output_window.see("end")


# display default values in the GUI
# def parameters(default_values, default_fit, default_constantF):
#     if not default_values == 0:
#         downsample_value.set(default_values['Downsampling rate'])
#         Filter_degree.set(default_values['Butterworth filter degree'])
#         Filter_cut_off.set(default_values['Cut-off frequency'])
#         Force_Min.set(default_values['Force threshold, pN'])
#         Z_score_force.set(default_values['Z-score force'])
#         Z_score_distance.set(default_values['Z-score distance'])
#         augment_factor_value.set(2)

#         step_d_variable.set(str(default_values['Step d']))
#         window_size_variable.set(str(default_values['Moving median window size']))
#         STD_difference_variable.set(str(default_values['STD difference threshold']))
#         Frequency_variable.set(str(default_values['Data frequency, Hz']))

#     dsLp_variable.set(str(default_fit['Persistance-Length ds, nm']))
#     dsLp_up_variable.set(str(default_fit['Persistance-Length ds, upper bound, nm']))
#     dsLp_low_variable.set(str(default_fit['Persistance-Length ds, lower bound, nm']))
#     ssLp_variable.set(str(default_fit['Persistance-Length ss, nm']))
#     dsLc_variable.set(str(default_fit['Contour-Length ds, nm']))
#     ssLc_variable.set(str(default_fit['Contour-Length ss, nm']))
#     ssLc_up_variable.set(str(default_fit['Contour-Length ss, upper bound, nm']))
#     stiff_ds_variable.set(str(default_fit['Stiffness ds, pN']))
#     stiff_ds_up_variable.set(str(default_fit['Stiffness ds, upper bound, pN']))
#     stiff_ds_low_variable.set(str(default_fit['Stiffness ds, lower bound, pN']))
#     stiff_ss_variable.set(str(default_fit['Stiffness ss, pN']))
#     stiff_ss_up_variable.set(str(default_fit['Stiffness ss, upper bound, pN']))
#     stiff_ss_low_variable.set(str(default_fit['Stiffness ss, lower bound, pN']))
#     f_off_variable.set(str(default_fit['Force offset, pN']))
#     f_off_up_variable.set(str(default_fit['Force offset, upper bound, pN']))
#     f_off_low_variable.set(str(default_fit['Force offset, lower bound, pN']))
#     d_off_variable.set(str(default_fit['Distance offset, nm']))
#     d_off_up_variable.set(str(default_fit['Distance offset, upper bound, nm']))
#     d_off_low_variable.set(str(default_fit['Distance offset, lower bound, nm']))

#     x_min.delete(0, "end")
#     x_min.insert("end", str(default_constantF['x min']))
#     x_max.delete(0, "end")
#     x_max.insert("end", str(default_constantF['x max']))
#     y_min.delete(0, "end")
#     y_min.insert("end", str(default_constantF['y min']))
#     y_max.delete(0, "end")
#     y_max.insert("end", str(default_constantF['y max']))
#     number_gauss.delete(0, "end")
#     number_gauss.insert("end", str(default_constantF['Number gauss']))
#     mean_gauss.delete(0, "end")
#     mean_gauss.insert("end", str(default_constantF['Mean']))
#     STD_gauss.delete(0, "end")
#     STD_gauss.insert("end", str(default_constantF['STD']))
#     amplitude_gauss.delete(0, "end")
#     amplitude_gauss.insert("end", str(default_constantF['Amplitude']))

def load_parameters():
    import_file_path = tk.filedialog.askopenfilename()  
    # Load the parameters from the text file
    with open(import_file_path, 'r') as file:
        lines = file.readlines()

    # Initialize dictionaries to hold the parameters
    default_values = {}
    default_fit = {}

    # Parse the parameters
    data_processing_start = lines.index("Data processing:\n")
    fitting_parameters_start = lines.index("Fitting parameters:\n")

    # Load Data processing parameters
    data_processing_lines = lines[data_processing_start + 1: fitting_parameters_start]
    default_values = json.loads(''.join(data_processing_lines).strip())

    # Load Fitting parameters
    fitting_parameters_lines = lines[fitting_parameters_start + 1:]
    default_fit = json.loads(''.join(fitting_parameters_lines).strip())

    # Set the GUI variables based on the loaded parameters
    downsample_value.set(default_values['downsample_value'])
    Filter_degree.set(default_values['filter_degree'])
    Filter_cut_off.set(default_values['filter_cut_off'])
    Force_Min.set(default_values['F_min'])
    Z_score_force.set(default_values['z-score_f'])
    Z_score_distance.set(default_values['z-score_d'])
    augment_factor_value.set(int(default_values['augment_factor']))

    step_d_variable.set(str(default_values['step_d']))
    window_size_variable.set(str(default_values['window_size']))
    STD_difference_variable.set(str(default_values['STD_diff']))
    Frequency_variable.set(str(default_values['data_frequency']))

    dsLp_variable.set(str(default_fit['lp_ds']))
    dsLp_up_variable.set(str(default_fit['lp_ds_up']))
    dsLp_low_variable.set(str(default_fit['lp_ds_low']))
    ssLp_variable.set(str(default_fit['lp_ss']))
    dsLc_variable.set(str(default_fit['lc_ds']))
    ssLc_variable.set(str(default_fit['lc_ss']))
    ssLc_up_variable.set(str(default_fit['lc_ss_up']))
    stiff_ds_variable.set(str(default_fit['ds_stiff']))
    stiff_ds_up_variable.set(str(default_fit['ds_stiff_up']))
    stiff_ds_low_variable.set(str(default_fit['ds_stiff_low']))
    stiff_ss_variable.set(str(default_fit['ss_stiff']))
    stiff_ss_up_variable.set(str(default_fit['ss_stiff_up']))
    stiff_ss_low_variable.set(str(default_fit['ss_stiff_low']))
    f_off_variable.set(str(default_fit['offset_f']))
    f_off_up_variable.set(str(default_fit['offset_f_up']))
    f_off_low_variable.set(str(default_fit['offset_f_low']))
    d_off_variable.set(str(default_fit['offset_d']))
    d_off_up_variable.set(str(default_fit['offset_d_up']))
    d_off_low_variable.set(str(default_fit['offset_d_low']))







# get all settings from the user input before start of the analysis
# def check_settings():
#     input_settings = {
#         'downsample_value': int(downsample_value2.get()),
#         'filter_degree': int(Filter_degree2.get()),
#         'filter_cut_off': float(Filter_cut_off2.get()),
#         'F_min': float(Force_Min2.get()),
#         'step_d': int(step_d_value.get()),
#         'z-score_f': float(Z_score_force2.get()),
#         'z-score_d': float(Z_score_distance2.get()),
#         'window_size': int(window_size_value.get()),
#         'data_frequency': float(Frequency_value.get()),
#         'STD_diff': float(STD_difference_value.get()),
#         'augment_factor': augment_factor_value.get()
#     }

#     input_format = {
#         'HF': check_box_HF.get(),
#         'LF': check_box_LF.get(),
#         'CSV': check_box_CSV.get(),
#         'Augment': check_box_augment.get(),
#         'Trap': check_box_Trap1.get(),
#         'length_measure': check_box_um.get(),
#         'MultiH5': check_box_multiH5.get(),
#         'preprocess': check_box_preprocess.get(),
#         'reverse_fitting':check_box_reverse_fitting.get()
#     }

#     export_data = {
#         'export_SMOOTH': check_box_smooth_data.get(),
#         'export_PLOT': check_box_plot.get(),
#         'export_STEPS': check_box_steps.get(),
#         'export_TOTAL': check_box_total_results.get(),
#         'export_FIT': check_box_fitting.get()
#     }

#     input_fitting = {
#         'WLC+WLC': int(check_box_WLC.get()),
#         'WLC+FJC': int(check_box_FJC.get()),
#         'lp_ds': float(dsLp.get()),
#         'lp_ds_up': float(dsLp_up.get()),
#         'lp_ds_low': float(dsLp_low.get()),
#         'lc_ds': float(dsLc.get()),
#         'lp_ss': float(ssLp.get()),
#         'lc_ss': float(ssLc.get()),
#         'lc_ss_up': float(ssLc_up.get()),
#         'ds_stiff': float(stiff_ds.get()),
#         'ds_stiff_up': float(stiff_ds_up.get()),
#         'ds_stiff_low': float(stiff_ds_low.get()),
#         'ss_stiff': float(stiff_ss.get()),
#         'ss_stiff_up': float(stiff_ss_up.get()),
#         'ss_stiff_low': float(stiff_ss_low.get()),
#         'offset_f': float(f_off.get()),
#         'offset_f_up': float(f_off_up.get()),
#         'offset_f_low': float(f_off_low.get()),
#         'offset_d': float(d_off.get()),
#         'offset_d_up': float(d_off_up.get()),
#         'offset_d_low': float(d_off_low.get())
#     }

#     input_constantF = {
#         'x min': int(x_min.get()),
#         'x max': int(x_max.get()),
#         'y min': int(y_min.get()),
#         'y max': int(y_max.get()),
#         'Number gauss': int(number_gauss.get()),
#         'Mean': mean_gauss.get(),
#         'STD': STD_gauss.get(),
#         'Amplitude': amplitude_gauss.get()
#     }

#     return input_settings, input_format, export_data, input_fitting, input_constantF


# export parameters used for the analysis in a txt file
def export_settings(analysis_path, timestamp, input_1, input_2):
    with open(str(analysis_path + '/parameters_' + timestamp + '.txt'), 'w') as config_used:
        config_used.write('Data processing:\n')
        config_used.write(json.dumps(input_1, indent=4, sort_keys=False))
        config_used.write('\n\n')
        config_used.write('Fitting parameters:\n')
        config_used.write(json.dumps(input_2, indent=4, sort_keys=False))


# Looks for output of the subprocess
def refresh():
    global new_image
    while output_q.empty() is False:
        output = output_q.get()
        output_window.insert("end", "\n" + output + "\n")
        output_window.see("end")
    try:
        images = str(analysis_folder + "/*plot*.png")
        list_images = glob.glob(images)
        img = Image.open(list_images[-1].replace("\\", "\\\\"))
        resized = img.resize((1000, 650), Image.ANTIALIAS)
        new_image = ImageTk.PhotoImage(resized)
        figure_frame.create_image((0, 0), image=new_image, anchor="nw")
    except:
        pass


def readme():
    with open("FRIED_POTATO_readme.txt", "r") as f:
        help_text = f.read()
        help_window = tk.Toplevel(root)
        help_window.title("Readme")
        text = tk.Text(help_window, height=25, width=200)
        text.grid(row=0, column=0, sticky="nw")
        text.insert("end", help_text)


# display a single file (tab2)
# def get_single_file():
#     # Specify the allowed file types
#     file_types = [
#         ("Text files", "*.txt"),
#         ("CSV files", "*.csv"),
#         ("h5 files", "*.h5")
#     ]

#     import_file_path = tk.filedialog.askopenfilename(file_types=file_types)

#     input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()

#     input_format['preprocess'] = 0
#     FD_raw, FD_raw_um, Frequency_value, filename = read_in_data(0, [import_file_path], input_settings, input_format)
#     input_format['preprocess'] = 1
#     FD, FD_um, Frequency_value, filename = read_in_data(0, [import_file_path], input_settings, input_format)
#     display_RAW_FD(FD[:, 0], FD[:, 1], FD_raw[:, 0], FD_raw[:, 1], filename)

#     if format == 'csv':
#         if not check_box_CSV.get() == 1:
#             check_box_CSV.set(value=1)
#             select_box(check_box_CSV, check_box_HF, check_box_LF)
#             parameters(default_values_CSV, default_values_FIT, default_values_constantF)
#         else:
#             pass


# create the plot for tab2
def display_RAW_FD(processed_F, processed_D, raw_F, raw_D, filename):
    global figure_raw

    try:
        figure_raw.get_tk_widget().destroy()
    except:
        pass

    single_fd = Figure(figsize=(10, 6), dpi=100)
    subplot1 = single_fd.add_subplot(111)

    legend_elements = [
        Line2D([0], [0], color='C0', lw=4),
        Line2D([0], [0], color='C1', lw=4)
    ]

    subplot1.set_title(str(filename))
    subplot1.set_xlabel("Distance (nm)")
    subplot1.set_ylabel("Force (pN)")
    subplot1.scatter(raw_D, raw_F, alpha=0.8, color='C0', s=0.1, zorder=0)
    subplot1.scatter(processed_D, processed_F, marker='.', s=0.1, linewidths=None, alpha=1, color='C1', zorder=1)
    subplot1.legend(legend_elements, ['Downsampled FD-Data', 'Filtered FD-Data'])

    figure_raw = FigureCanvasTkAgg(single_fd, figure_frame2)
    figure_raw.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab2)


def start_constantF():
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    Force_Distance, Force_Distance_um, frequency, filename, analysis_path, timestamp = get_constantF(input_settings, input_format, input_constantF)
    fig_constantF, hist_D, filteredDistance_ready = display_constantF(Force_Distance, Force_Distance_um, frequency, input_settings, input_constantF)
    os.mkdir(analysis_path)
    export_settings(analysis_path, timestamp, input_settings, input_constantF)
    fig_constantF = fit_constantF(hist_D, Force_Distance, filteredDistance_ready, frequency, input_settings, input_constantF, filename, timestamp)
    fig_constantF_tk = FigureCanvasTkAgg(fig_constantF, figure_frame_tab4)
    fig_constantF_tk.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab4)


def show_constantF():
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    Force_Distance, Force_Distance_um, frequency, filename, analysis_path, timestamp = get_constantF(input_settings, input_format, input_constantF)
    fig_constantF, hist_D, filteredDistance_ready = display_constantF(Force_Distance, Force_Distance_um, frequency, input_settings, input_constantF)
    fig_constantF_tk = FigureCanvasTkAgg(fig_constantF, figure_frame_tab4)
    fig_constantF_tk.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab4)





################ TOMATO ###############################
from FRIED_POTATO_TOMATO import plot_TOMATO


############# define the functions for TOMATO ##################
def open_folder():
    global filename_TOMATO
    global Force_Distance_TOMATO
    global der_arr_TOMATO
    global TOMATO_fig1
    global Files
    global FD_number
    # check user input
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()

    # ask wich directory should be analysed
    folder = tk.filedialog.askdirectory()
    root.title('POTATO -- ' + str(folder))

    # decide which input format was choosen
    if input_format['CSV'] == 1:
        folder_path = str(folder + "/*.csv")
    else:
        folder_path = str(folder + "/*.h5")

    Files = glob.glob(folder_path)

    FD_number = 0
    Force_Distance_TOMATO, Force_Distance_um_TOMATO, Frequency_value, filename_TOMATO = read_in_data(FD_number, Files, input_settings, input_format)
    der_arr_TOMATO = create_derivative(input_settings, Frequency_value, Force_Distance_TOMATO[:, 0], Force_Distance_TOMATO[:, 1], 0)

    entryText_filename.set(filename_TOMATO)

    try:
        TOMATO_fig1.get_tk_widget().destroy()
    except:
        pass

    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_figure_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_figure_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def change_FD(direction):
    global TOMATO_fig1
    global filename_TOMATO
    global FD_number
    global Force_Distance_TOMATO
    global orientation

    FD_number = FD_number + direction
    if FD_number == len(Files):
        FD_number = FD_number - len(Files)
    if FD_number < 0:
        FD_number = FD_number + len(Files)

    delete_all_steps()
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    Force_Distance_TOMATO, Force_Distance_um_TOMATO, Frequency_value, filename_TOMATO = read_in_data(FD_number, Files, input_settings, input_format)

    orientation = 'forward'
    if Force_Distance_TOMATO[0, 1] > Force_Distance_TOMATO[-1, 1]:  # reverse
        Force_Distance_TOMATO = np.flipud(Force_Distance_TOMATO)
        Force_Distance_um_TOMATO = np.flipud(Force_Distance_um_TOMATO)
        orientation = 'reverse'

    entryText_filename.set(filename_TOMATO)

    parameters(0, default_values_FIT, default_values_constantF)

    TOMATO_fig1.get_tk_widget().destroy()

    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_figure_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_figure_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


# key binding wrapper functions
def previous_FD_key(event):
    change_FD(-1)


def next_FD_key(event):
    change_FD(+1)


def save_step_key(event):
    save_step()


def start_analysis_key(event):
    analyze_steps()


def start_click_key(event):
    start_click()


def end_click_key(event):
    end_click()


def start_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=1: onclick_start_end(event, arg))


def end_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=0: onclick_start_end(event, arg))


def onclick_start_end(event, pos):
    global cid

    PD_position, F_position = float(event.xdata), float(event.ydata)

    if pos == 1:
        entryText_startF.set(round(F_position, 1))
        entryText_startD.set(round(PD_position, 1))
    elif pos == 0:
        entryText_endF.set(round(F_position, 1))
        entryText_endD.set(round(PD_position, 1))
    TOMATO_fig1.mpl_disconnect(cid)


def save_step():
    global step_number
    try:
        tree_steps.insert('', 'end', values=(step_number, entryText_startF.get(), entryText_endF.get(), entryText_startD.get(), entryText_endD.get()))
        step_number += 1
    except:
        print('Please make sure step start and step end are selected!')


def analyze_steps():
    global TOMATO_fig1
    global subplot1
    global tree_results

    # get input settings
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    timestamp = time.strftime("%Y%m%d-%H%M%S")

    # write step list into pandas dataframe
    row_list = []
    columns = ('Step number', 'F start', 'F end', 'Step start', 'Step end')
    for row in tree_steps.get_children():
        row_list.append(tree_steps.item(row)["values"])
    treeview_df = pd.DataFrame(row_list, columns=columns)

    # iterate through dataframe and fit each part of the curve
    TOMATO_fig1.get_tk_widget().destroy()

    figure1 = plot_TOMATO(Force_Distance_TOMATO)
    diff_colors = ['b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm']
    subplot1 = figure1.add_subplot(111)
    subplot1.plot(Force_Distance_TOMATO[:, 1], Force_Distance_TOMATO[:, 0], color='gray')
    distance = np.arange(min(Force_Distance_TOMATO[:, 1]), max(Force_Distance_TOMATO[:, 1]) + 50, 2)

    export_fit = []
    fit = []
    start_force_ss = []
    start_distance_ss = []
    integral_ss_fit_start = []
    integral_ss_fit_end = []

    for i in treeview_df.index:
        # part before first step is fitted with a single WLC model (ds part)
        if treeview_df['Step number'][i] == 1:
            j = treeview_df['Step number'][i]
            ds_fit_dict_TOMATO, TOMATO_area_ds, real_step_start = fitting_ds(filename_TOMATO, input_settings, export_data, input_fitting, float(treeview_df['Step start'][i]), Force_Distance_TOMATO, der_arr_TOMATO, [], 1)
            ds_fit_region_end = real_step_start

            dsLp_variable.set(ds_fit_dict_TOMATO['Lp_ds'])
            f_off_variable.set(ds_fit_dict_TOMATO['f_offset'])
            d_off_variable.set(ds_fit_dict_TOMATO["d_offset"])
            dsLc_variable.set(ds_fit_dict_TOMATO['Lc_ds'])
            stiff_ds_variable.set(ds_fit_dict_TOMATO['St_ds'])

            tree_results.insert("", "end", iid='{}no step'.format(timestamp), values=(
                entryText_filename.get(),
                i,
                '',
                '',
                '',
                '',
                '',
                '',
                dsLc_variable.get(),
                dsLp_variable.get(),
                stiff_ds_variable.get(),
                '',
                '',
                '',
                f_off_variable.get(),
                d_off_variable.get(),
                '',
                ''
            )
            )

            export_fit.append(ds_fit_dict_TOMATO)

            F_ds_model = ds_fit_dict_TOMATO['model_ds'](distance, ds_fit_dict_TOMATO['fit_model'].params)
            # plot the marked ds region and fits
            subplot1.plot(Force_Distance_TOMATO[:, 1][:real_step_start], Force_Distance_TOMATO[:, 0][:real_step_start], color=diff_colors[i])
            subplot1.plot(distance, F_ds_model, marker=None, linestyle='dashed', linewidth=1, color="black")

        # fit the other ss parts
        elif treeview_df['Step number'][i] > 1:
            j = treeview_df['Step number'][i]
            fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict_TOMATO, area_ss_fit_start, area_ss_fit_end = fitting_ss(filename_TOMATO, input_settings, export_data, input_fitting, float(treeview_df['Step end'][i - 1]), float(treeview_df['Step start'][i]), Force_Distance_TOMATO, 1, 1, der_arr_TOMATO, [], 1)

            fit.append(fit_ss)
            start_force_ss.append(f_fitting_region_ss)
            start_distance_ss.append(d_fitting_region_ss)
            export_fit.append(ss_fit_dict_TOMATO)
            integral_ss_fit_start.append(area_ss_fit_start)
            integral_ss_fit_end.append(area_ss_fit_end)

            ssLp_variable.set(ss_fit_dict_TOMATO['Lp_ss'])
            f_off_variable.set(ss_fit_dict_TOMATO['f_offset'])
            d_off_variable.set(ss_fit_dict_TOMATO["d_offset"])
            ssLc_variable.set(ss_fit_dict_TOMATO['Lc_ss'])
            stiff_ss_variable.set(ss_fit_dict_TOMATO['St_ss'])

            tree_results.insert("", "end", iid='{}step{}'.format(timestamp, j-1), values=(
                entryText_filename.get(),
                i,
                Force_Distance_TOMATO[real_step_start, 0],
                f_fitting_region_ss[0],
                (f_fitting_region_ss[0] + Force_Distance_TOMATO[real_step_start, 0]) / 2,
                Force_Distance_TOMATO[real_step_start, 1],
                d_fitting_region_ss[0],
                d_fitting_region_ss[0] - Force_Distance_TOMATO[real_step_start, 1],
                '',
                '',
                '',
                ssLc_variable.get(),
                ssLp_variable.get(),
                stiff_ss_variable.get(),
                f_off_variable.get(),
                d_off_variable.get(),
                '',
                ''
            )
            )

            real_step_start = np.where(Force_Distance_TOMATO[:, 0] == f_fitting_region_ss[-1])
            real_step_start = real_step_start[0][0]

            # plot the marked regions and fits
            # model data
            F_ss_model = ss_fit_dict_TOMATO['model_ss'](distance, fit_ss.params)

            # plot the marked ss region and fits
            subplot1.plot(d_fitting_region_ss[:], f_fitting_region_ss, color=diff_colors[i])
            subplot1.plot(distance, F_ss_model, marker=None, linewidth=1, linestyle='dashed', color="black")

    # fit the last part of the curve
    fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict_TOMATO, area_ss_fit_start, area_ss_fit_end = fitting_ss(
        filename_TOMATO,
        input_settings,
        export_data,
        input_fitting,
        float(treeview_df['Step end'][len(treeview_df) - 1]),
        max(Force_Distance_TOMATO[:, 1]),
        Force_Distance_TOMATO,
        1,
        1,
        der_arr_TOMATO,
        [],
        1
    )

    fit.append(fit_ss)
    start_force_ss.append(f_fitting_region_ss)
    start_distance_ss.append(d_fitting_region_ss)
    export_fit.append(ss_fit_dict_TOMATO)
    integral_ss_fit_start.append(area_ss_fit_start)
    integral_ss_fit_end.append(area_ss_fit_end)

    ssLp_variable.set(ss_fit_dict_TOMATO['Lp_ss'])
    f_off_variable.set(ss_fit_dict_TOMATO['f_offset'])
    d_off_variable.set(ss_fit_dict_TOMATO["d_offset"])
    ssLc_variable.set(ss_fit_dict_TOMATO['Lc_ss'])
    stiff_ss_variable.set(ss_fit_dict_TOMATO['St_ss'])

    tree_results.insert("", "end", iid='{}step{}'.format(timestamp, j), values=(
        entryText_filename.get(),
        j,
        Force_Distance_TOMATO[:, 0][real_step_start],
        f_fitting_region_ss[0],
        (f_fitting_region_ss[0] + Force_Distance_TOMATO[:, 0][real_step_start]) / 2,
        Force_Distance_TOMATO[:, 1][real_step_start],
        d_fitting_region_ss[0],
        d_fitting_region_ss[0] - Force_Distance_TOMATO[:, 1][real_step_start],
        '',
        '',
        '',
        ssLc_variable.get(),
        ssLp_variable.get(),
        stiff_ss_variable.get(),
        f_off_variable.get(),
        d_off_variable.get(),
        '',
        ''
    )
    )

    work_first_step, kT_1 = calc_integral(
        TOMATO_area_ds,
        integral_ss_fit_start[0],
        Force_Distance_TOMATO[ds_fit_region_end, 1],
        start_distance_ss[0][0],
        Force_Distance_TOMATO[ds_fit_region_end, 0],
        start_force_ss[0][0]
    )

    tree_results.set('{}step1'.format(timestamp), column='Work [pN*nm]', value=work_first_step)
    tree_results.set('{}step1'.format(timestamp), column='Work [kT]', value=kT_1)

    if j > 1:
        for n in range(1, j):
            print(start_distance_ss[n - 1][-1])
            print(start_distance_ss[n][0])
            work_step_n, kT_n = calc_integral(
                integral_ss_fit_end[n - 1],
                integral_ss_fit_start[n],
                start_distance_ss[n - 1][-1],
                start_distance_ss[n][0],
                start_force_ss[n - 1][-1],
                start_force_ss[n][0]
            )
            print('WORK', work_step_n)
            tree_results.set('{}step{}'.format(timestamp, n+1), column='Work [pN*nm]', value=work_step_n)
            tree_results.set('{}step{}'.format(timestamp, n+1), column='Work [kT]', value=kT_n)

    # plot the marked regions and fits
    # model data
    F_ss_model = ss_fit_dict_TOMATO['model_ss'](distance, fit_ss.params)

    # plot the marked ss region and fits
    subplot1.plot(d_fitting_region_ss[:], f_fitting_region_ss, color=diff_colors[j + 1])
    subplot1.plot(distance, F_ss_model, marker=None, linewidth=1, linestyle='dashed', color="black")

    subplot1.set_ylim([min(Force_Distance_TOMATO[:, 0]), max(Force_Distance_TOMATO[:, 0])])
    subplot1.set_xlim([min(Force_Distance_TOMATO[:, 1]) - 10, max(Force_Distance_TOMATO[:, 1]) + 10])
    subplot1.tick_params('both', direction='in')

    TOMATO_fig1 = FigureCanvasTkAgg(figure1, TOMATO_figure_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0)

    toolbarFrame = tk.Frame(master=TOMATO_figure_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def delete_step():
    global step_number
    list_items = tree_steps.get_children("")
    tree_steps.delete(list_items[-1])
    step_number -= 1


def delete_all_steps():
    global step_number

    tree_steps.delete(*tree_steps.get_children())
    step_number = 1


def delete_result(event):
    selected_items = tree_results.selection()
    for selected_item in selected_items:
        tree_results.delete(selected_item)


def clear_table():
    global tree_results
    list_items = tree_results.get_children("")

    for item in list_items:
        tree_results.delete(item)


def clear_table_last():
    global tree_results
    list_items = tree_results.get_children("")

    tree_results.delete(list_items[-1])


def export_table():
    global tree_results
    global name
    global Fit_results
    ''' exporting the table results '''
    results = []
    for child in tree_results.get_children():
        results.append(tree_results.item(child)['values'])

    Fit_results = pd.DataFrame(results, columns=[
        'Filename',
        'step number',
        'Force step start [pN]',
        'Force step end [pN]',
        'mean force [pN]',
        'extension step start [nm]',
        'extension step end [nm]',
        'Step length [nm]',
        'ds contour length',
        'ds persistance Length',
        'ds stiffness (K0) [pN]',
        'ss contour Length',
        'ss persistance Length',
        'ss stiffness (K0) [pN]',
        'Force offset',
        'Distance offset',
        'Work [pN*nm]',
        'Work [kT]'
        ])

    name = tk.filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    Fit_results.to_csv(name.name, index=False, header=True)


def tab_bind(event=None):
    if tabControl.index(tabControl.select()) == 4:
        root.bind("<Right>", next_FD_key)
        root.bind("<Left>", previous_FD_key)
        root.bind("<s>", start_click_key)
        root.bind("<e>", end_click_key)
        root.bind("<Control-s>", save_step_key)
        root.bind("<Control-f>", start_analysis_key)
        root.bind("<Delete>", delete_result)
    else:
        root.unbind("<Right>")
        root.unbind("<Left>")
        root.unbind("<s>")
        root.unbind("<e>")
        root.unbind("<Control-s>")
        root.unbind("<Control-f>")
        root.unbind("<Delete>")

############## TOMATO functions end ###################

"""class defining the tkinter GUI"""
class FriedPotatoGUI:
    def __init__(self):
        # create the main window
        mp.freeze_support()
        self.root = tk.Tk()
        self.root.minsize(1300, 800)
        self.root.iconbitmap('FRIED_POTATO.ico')
        self.root.title("FRIED POTATO -- Force-Ramp Improved EDition of Practical Optical Tweezers Analysis TOol")
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # create a queue for the text output of the subprocesses
        self.output_q = mp.Queue()

        # Initialize checkboxes and checkbox_vars as dictionaries to store checkbox objects and their state variables
        self.checkboxes = {}
        self.checkbox_vars = {}

        # Initialize text input fields and their variables as dictionaries to store text input objects and their state variables
        self.text_input = {}
        self.text_input_vars = {}
    
        # setup the drop-down menus and the different tabs
        self.setup_menus()
        self.setup_tabs()

        # create a text output window that is shared between all tabs
        self.display_text_output()

        # basic layout for all tabs, but the settings
        self.tab_frames = {}
        self.divide_tab(self.tab_main, 'main')
        self.divide_tab(self.tab_display_curve, 'display_curve')
        self.divide_tab(self.tab_constant_force, 'constant_force')
        self.divide_tab(self.tab_TOMATO, 'TOMATO')
        self.root.bind('<<NotebookTabChanged>>', self.tab_bind)

        # specific layout for the main tab
        self.layout_tab_main()
        self.layout_tab_display_curve()
        self.layout_tab_settings()
        self.layout_tab_constant_force()
        self.layout_tab_TOMATO()

        # put default values into the widgets
        self.parameters(default_values_HF, default_values_FIT, default_values_constantF)


    def setup_menus(self):
        drop_down_menu = tk.Menu(self.root)
        self.root.config(menu=drop_down_menu)

        # File menu
        file_menu = tk.Menu(drop_down_menu, tearoff=0)
        drop_down_menu.add_cascade(label='File', menu=file_menu)
        file_menu.add_command(label='Analyse folder (FD curves)', command=self.start_analysis)
        file_menu.add_command(label='Display single FD curve', command=lambda: self.get_single_file())
        file_menu.add_command(label='Show h5 file structure', command=lambda: self.show_h5())
        file_menu.add_command(label="Load Parameters", command=lambda: self.load_parameters())
        file_menu.add_separator()
        file_menu.add_command(label='Display constant force', command=self.show_constantF)
        file_menu.add_command(label='Fit constant force', command=self.start_constantF)

        # Settings menu
        settings_menu = tk.Menu(drop_down_menu, tearoff=0)
        drop_down_menu.add_cascade(label='Settings', menu=settings_menu)
        settings_menu.add_command(label='Set advanced settings', command=lambda: self.tabControl.select(self.tab3))

        # Help menu
        help_menu = tk.Menu(drop_down_menu, tearoff=0)
        drop_down_menu.add_cascade(label='Help', menu=help_menu)
        help_menu.add_command(label='Readme', command=self.readme)


    def display_text_output(self):
        # Create a scrollable output window
        self.output_window = tk.Text(self.root, height=10, width=115)
        self.output_window.grid(row=0, column=0, sticky="nsew")
    
        # Create a scrollbar
        scrollbar = tk.Scrollbar(self.root, command=self.output_window.yview)
        scrollbar.grid(row=0, column=1)
    
        self.output_window.config(yscrollcommand=scrollbar.set)
    
        self.output_window.insert(
            "end",
            "Welcome to FRIED POTATO! \n"
            "This tool can analyze force-distance curves obtained from optical tweezers experiments. \n"
            "This can be done either in 'bulk' mode, where all files in a folder are analysed, or in 'manual' TOMATO mode. \n"
            "Parameters should be adjusted prior to analysis.\n"
        )
        self.output_window.config(state="disabled")  # Make the window read-only

        # add a refresh button to the output window
        self.create_button('Refresh', self.root, self.refresh, 0, 2, 5, 2, 3, 6, cursor="exchange")


    def setup_tabs(self):
        self.style = ttk.Style()

        mybackground = "#b5b8b5"
        myforeground = "#8c888c"

        self.style.theme_create("fried", parent="alt", settings={
            "TNotebook": {"configure": {"tabmargins": [2, 5, 2, 0] }},
            "TNotebook.Tab": {
                "configure": {"padding": [5, 1], "background": mybackground},
                "map":       {"background": [("selected", myforeground)],
                "expand":    [("selected", [1, 1, 1, 0])] } 
            }
        })

        self.style.theme_use("fried")

        self.tabControl = ttk.Notebook(self.root)
        self.tabControl.grid(row=1, column=0, padx=2, pady=2, columnspan=3, sticky='nsew')

        self.root.grid_rowconfigure(0, weight=1)  # Make the row containing the tabControl expandable
        self.root.grid_columnconfigure(0, weight=1)  # Make the column containing the tabControl expandable

        self.tab_main = ttk.Frame(self.tabControl, width=800, height=600)
        self.tab_display_curve = ttk.Frame(self.tabControl, width=800, height=600)
        self.tab_settings = ttk.Frame(self.tabControl, width=800, height=600)
        self.tab_constant_force = ttk.Frame(self.tabControl, width=800, height=600)
        self.tab_TOMATO = ttk.Frame(self.tabControl, width=800, height=600)

        self.tabControl.add(self.tab_main, text="Folder Analysis")
        self.tabControl.add(self.tab_display_curve, text="Show Single File")
        self.tabControl.add(self.tab_settings, text="Advanced Settings")
        self.tabControl.add(self.tab_constant_force, text="Constant Force Analysis")
        self.tabControl.add(self.tab_TOMATO, text="Manual Analysis - TOMATO")


    def divide_tab(self, tab, tab_name):
        """ divide the tab into frames """
        output_frame = tk.Canvas(tab, height=550, width=800, borderwidth=1)

        output_frame.grid(row=0, column=0, sticky="nsew")
        tab.grid_columnconfigure(0, weight=1)  # Make the column expandable
        tab.grid_rowconfigure(0, weight=1)  # Make the row expandable

        # Create a separator
        separator = ttk.Separator(tab, orient='vertical')
        separator.grid(row=0, column=1, sticky='ns')

        input_frame = tk.Frame(tab, width=200, height=500)
        input_frame.grid(row=0, column=2, sticky="nsew")
        tab.grid_columnconfigure(2, weight=1)  # Make the column expandable



        # Store the frames in the dictionary using the tab_name as the key
        self.tab_frames[tab_name] = {"output_frame": output_frame, "input_frame": input_frame}


    def layout_tab_main(self):
        output_frame = self.tab_frames["main"]["output_frame"]
        input_frame = self.tab_frames["main"]["input_frame"]

        # input parameters - check boxes
        check_box_frame = tk.Frame(input_frame)
        check_box_frame.grid(row=0, column=0)

        # check boxes to specify the data type
        cluster_data_type = tk.Label(check_box_frame, text='DATA TYPE', font='Helvetica 9 bold')
        cluster_data_type.grid(row=0, column=0, padx=2, pady=(2, 2), sticky='W')

        self.create_check_box("HF", check_box_frame, "High Frequency (Piezo Distance)", 1, lambda: self.select_box("HF", "LF", "CSV"), 1, 0, 0)
        self.create_check_box("LF", check_box_frame, "Low Frequency", 0, lambda: self.select_box("LF", "HF", "CSV"), 2, 0, 0)
        self.create_check_box("CSV", check_box_frame, "CSV (F(pN) | d)", 0, lambda: self.select_box('CSV', 'HF', 'LF'), 3, 0, 0)

        # check boxes to specify the data handling
        cluster_handling = tk.Label(check_box_frame, text='HANDLING', font='Helvetica 9 bold')
        cluster_handling.grid(row=4, column=0, padx=2, pady=(10, 2), sticky='W')

        self.create_check_box("multi", check_box_frame, "MultiH5", 0, None, 5, 0, 0)
        self.create_check_box("rev", check_box_frame, "Reverse Fitting", 0, None, 6, 0, 0)

        # check boxes to specify the data location
        cluster_data_location = tk.Label(check_box_frame, text='DATA LOCATION', font='Helvetica 9 bold')
        cluster_data_location.grid(row=0, column=1, padx=2, pady=(2, 2), sticky='W')

        self.create_check_box("1x", check_box_frame, "Trap 1x", 0, lambda: self.select_box("1x", "2x"), 1, 1, 0)
        self.create_check_box("2x", check_box_frame, "Trap 2x", 1, lambda: self.select_box("2x", "1x"), 2, 1, 0)

        # check boxes to specify the data units
        cluster_data_units = tk.Label(check_box_frame, text='DATA UNITS', font='Helvetica 9 bold')
        cluster_data_units.grid(row=4, column=1, padx=2, pady=(10, 2), sticky='W')

        self.create_check_box("um", check_box_frame, "µm input", 1, lambda: self.select_box("um", "nm"), 5, 1, 0)
        self.create_check_box("nm", check_box_frame, "nm input", 0, lambda: self.select_box("nm", "um"), 6, 1, 0)


        # input parameters - text
        parameter_frame = tk.Frame(input_frame)
        parameter_frame.grid(row=1, column=0, sticky='NE')

        cluster_preprocessing = tk.Label(parameter_frame, text='PREPROCESSING', font='Helvetica 9 bold')
        cluster_preprocessing.grid(row=0, column=0, padx=2, pady=(20, 2))
        self.create_check_box("prepro", parameter_frame, "Preprocessing", 1, None, 0, 1, (20, 2))
        self.create_text_input("downsample", parameter_frame, 'Downsampling rate', 1, 0, 2, 2)
        self.create_text_input("filter_degree", parameter_frame, 'Butterworth filter degree', 2, 0, 2, 2)
        self.create_text_input("filter_cut_off", parameter_frame, 'Cut-off frequency', 3, 0, 2, 2)
        self.create_text_input("force_min", parameter_frame, 'Force threshold, pN', 4, 0, 2, 2)

        cluster_statistics = tk.Label(parameter_frame, text='STATISTICS', font='Helvetica 9 bold')
        cluster_statistics.grid(row=5, column=0, padx=2, pady=(20, 2))
        self.create_text_input("z_score_force", parameter_frame, 'Z-score force', 6, 0, 2, 2)
        self.create_text_input("z_score_distance", parameter_frame, 'Z-score distance', 7, 0, 2, 2)

        cluster_augment = tk.Label(parameter_frame, text='AUGMENTATION', font='Helvetica 9 bold')
        cluster_augment.grid(row=8, column=0, padx=2, pady=(20, 2))
        self.create_check_box("aug", parameter_frame, "Data Augmentation", 0, None, 8, 1, (20, 2))
        self.create_text_input("augment_factor", parameter_frame, 'Augmentation factor', 9, 0, 2, 2)

        self.create_button("Select Folder to Analyse!", parameter_frame, self.start_analysis, 11, 0, 2, 50, 2, 20, 2)


    def layout_tab_display_curve(self):
        output_frame = self.tab_frames["display_curve"]["output_frame"]
        input_frame = self.tab_frames["display_curve"]["input_frame"]

        open_file = tk.Button(
            input_frame,
            text='Open file to display',
            command=lambda: self.get_single_file(),
            bg='#df4c4c',
            activebackground='#eaa90d',
            font='Helvetica 11 bold',
            height=1,
            width=15
        )

        open_file.grid(row=0, column=0, pady=20, sticky='nsew')


    def layout_tab_settings(self):
        # Configure the grid to expand
        self.tab_settings.grid_columnconfigure(0, weight=1)
        self.tab_settings.grid_columnconfigure(2, weight=1)
        self.tab_settings.grid_columnconfigure(4, weight=1)
        self.tab_settings.grid_rowconfigure(0, weight=1)

        step_finding_settings = tk.Canvas(self.tab_settings, width=300, height=600)
        step_finding_settings.grid(row=0, column=0, sticky='NSEW')

        separator1 = ttk.Separator(self.tab_settings, orient='vertical')
        separator1.grid(row=0, column=1, sticky='NS')

        export_settings = tk.Canvas(self.tab_settings, width=300, height=600)
        export_settings.grid(row=0, column=2, sticky='NSEW')

        separator2 = ttk.Separator(self.tab_settings, orient='vertical')
        separator2.grid(row=0, column=3, sticky='NS')

        fitting_settings = tk.Canvas(self.tab_settings, width=500, height=600)
        fitting_settings.grid(row=0, column=4, sticky='NSEW')

        # """ step finding settings """
        cluster_preprocessing = tk.Label(step_finding_settings, text='PREPROCESSING', font='Helvetica 9 bold')
        cluster_preprocessing.grid(row=0, column=0, padx=2, pady=(20, 2))
        self.create_text_input("downsample", step_finding_settings, 'Downsampling rate', 1, 0, 2, 2)
        self.create_text_input("filter_degree", step_finding_settings, 'Butterworth filter degree', 2, 0, 2, 2)
        self.create_text_input("filter_cut_off", step_finding_settings, 'Cut-off frequency', 3, 0, 2, 2)
        self.create_text_input("force_min", step_finding_settings, 'Force threshold, pN', 4, 0, 2, 2)

        cluster_derivative = tk.Label(step_finding_settings, text="DERIVATIVE", font='Helvetica 9 bold')
        cluster_derivative.grid(row=5, column=0, padx=2, pady=(20, 2))
        self.create_text_input("step_d", step_finding_settings, 'Step d', 6, 0, 2, 2)
        self.create_text_input("frequency", step_finding_settings, 'Data frequency, Hz', 7, 0, 2, 2)

        cluster_statistics = tk.Label(step_finding_settings, text='STATISTICS', font='Helvetica 9 bold')
        cluster_statistics.grid(row=8, column=0, padx=2, pady=(20, 2))
        self.create_text_input("z_score_force", step_finding_settings, 'Z-score force', 9, 0, 2, 2)
        self.create_text_input("z_score_distance", step_finding_settings, 'Z-score distance', 10, 0, 2, 2)
        self.create_text_input("window_size", step_finding_settings, 'Moving median window size', 12, 0, 2, 2)
        self.create_text_input("std_difference", step_finding_settings, 'SD difference threshold', 13, 0, 2, 2)

        # """ Export settings """
        cluster_export = tk.Label(export_settings, text='EXPORT', font='Helvetica 9 bold')
        cluster_export.grid(row=0, column=0, padx=20, pady=20)
        self.create_check_box("smooth_data", export_settings, "Smoothed data", 1, None, 1, 0, 0)
        self.create_check_box("plot", export_settings, "Plot", 1, None, 2, 0, 0)
        self.create_check_box("steps", export_settings, "Steps found", 1, None, 3, 0, 0)
        self.create_check_box("total_results", export_settings, "Total results", 1, None, 4, 0, 0)
        self.create_check_box("fitting", export_settings, "Fitting", 1, None, 5, 0, 0)

        # """ Fitting parameters """
        # # Labels
        cluster_fitting = tk.Label(fitting_settings, text='FITTING', font='Helvetica 9 bold')
        cluster_fitting.grid(row=0, column=0, padx=20, pady=20)
        self.create_check_box("WLC", fitting_settings, "WLC+WLC", 1, lambda: self.select_box("WLC", "FJC"), 1, 0, 2)
        self.create_check_box("FJC", fitting_settings, "WLC+FJC", 0, lambda: self.select_box("FJC", "WLC"), 2, 0, 2)
        self.create_text_input("dsLp", fitting_settings, 'dsLp [nm]', 3, 0, 2, (40,2))
        self.create_text_input("dsLp_up", fitting_settings, 'dsLp upper bound [nm]', 4, 0, 2, 2)
        self.create_text_input("dsLp_low", fitting_settings, 'dsLp lower bound [nm]', 5, 0, 2, 2)
        self.create_text_input("dsLc", fitting_settings, 'dsLc [nm]', 6, 0, 2, (40,2))
        self.create_text_input("ssLp", fitting_settings, 'ssLp [nm]', 3, 2, (20,2), (40,2))
        self.create_text_input("ssLc", fitting_settings, 'ssLc [nm]', 6, 2, (20,2), (40,2))
        self.create_text_input("ssLc_up", fitting_settings, 'ssLc upper bound [nm]', 7, 2, (20,2), 2)
        self.create_text_input("stiffness_ds", fitting_settings, 'dsK0 [pN]', 8, 0, 2, (40,2))
        self.create_text_input("stiffness_ds_up", fitting_settings, 'dsK0 upper bound [pN]', 9, 0, 2, 2)
        self.create_text_input("stiffness_ds_low", fitting_settings, 'dsK0 lower bound [pN]', 10, 0, 2, 2)
        self.create_text_input("stiffness_ss", fitting_settings, 'ssK0 [pN]', 8, 2, (20,2), (40,2))
        self.create_text_input("stiffness_ss_up", fitting_settings, 'ssK0 upper bound [pN]', 9, 2, (20,2), 2)
        self.create_text_input("stiffness_ss_low", fitting_settings, 'ssK0 lower bound [pN]', 10, 2, (20,2), 2)
        self.create_text_input("f_offset", fitting_settings, 'Force offset [pN]', 11, 0, 2, (40,2))
        self.create_text_input("f_offset_up", fitting_settings, 'Force offset upper bound [pN]', 12, 0, 2, 2)
        self.create_text_input("f_offset_low", fitting_settings, 'Force offset lower bound [pN]', 13, 0, 2, 2)
        self.create_text_input("d_offset", fitting_settings, 'Distance offset [nm]', 11, 2, (20,2), (40,2))
        self.create_text_input("d_offset_up", fitting_settings, 'Distance offset upper bound [nm]', 12, 2, (20,2), 2)
        self.create_text_input("d_offset_low", fitting_settings, 'Distance offset lower bound [nm]', 13, 2, (20,2), 2)


    def layout_tab_constant_force(self):
        output_frame = self.tab_frames["constant_force"]["output_frame"]
        input_frame = self.tab_frames["constant_force"]["input_frame"]

        self.create_button("Fit Constant Force Data", input_frame, self.start_constantF, 0, 0, 20, 20, 2, 25, 2)
        self.create_button("Display Constant Force Data", input_frame, self.show_constantF, 1, 0, 20, 20, 2, 25, 2)

        # organize settings
        cluster_axes = tk.Label(input_frame, text='SET AXES', font='Helvetica 9 bold')
        cluster_axes.grid(row=3, column=0, padx=2, pady=(20, 2))
        self.create_text_input("x_min_constant_f", input_frame, 'x min', 4, 0, 2, 2)
        self.create_text_input("x_max_constant_f", input_frame, 'x max', 5, 0, 2, 2)

        self.create_text_input("y_min_constant_f", input_frame, 'y min', 4, 2, (20, 0), 2)
        self.create_text_input("y_max_constant_f", input_frame, 'y max', 5, 2, (20, 0), 2)

        cluster_expected_fit = tk.Label(input_frame, text='EXPECTED VALUES', font='Helvetica 9 bold')
        cluster_expected_fit.grid(row=6, column=0, padx=2, pady=(20, 2))
        self.create_text_input("number_gauss", input_frame, 'Number of expected gaussians', 7, 0, 2, 2)
        self.create_text_input("mean_gauss", input_frame, 'Expected mean of each gaussian', 8, 0, 2, 2)
        self.create_text_input("STD_gauss", input_frame, 'Expected SD of each gaussian', 9, 0, 2, 2)
        self.create_text_input("amplitude_gauss", input_frame, 'Expected amplitude of each gaussian', 10, 0, 2, 2)


    def layout_tab_TOMATO(self):
        output_frame = self.tab_frames["TOMATO"]["output_frame"]
        input_frame = self.tab_frames["TOMATO"]["input_frame"]

        # split output frame into two parts
        TOMATO_figure_frame = tk.Canvas(output_frame, width=800, height=590)
        TOMATO_figure_frame.grid(row=0, column=0, sticky='nsew')

        TOMATO_results_frame = tk.Canvas(output_frame, width=800, height=150)
        TOMATO_results_frame.grid(row=1, column=0, sticky='nsew')

        ### create buttons ###
        self.create_button("Choose folder", input_frame, self.start_TOMATO, 0, 0, 4, 4, 2, 16, 2)
        self.create_button("Save results table", input_frame, self.export_TOMATO, 1, 0, 4, 4, 2, 16, 1)
        self.create_button("Reset parameters", input_frame, self.reset_TOMATO, 1, 1, 4, 4, 2, 16, 1)
    
        # ### create entry widgets ###
        # inital shift in distance and force
        self.create_text_input("d_offset", input_frame, 'Distance offset [nm]', 2, 0, 2, 2)
        self.create_text_input("f_offset", input_frame, 'Force offset [pN]', 2, 2, 2, 2)

        # K0 for ds and ss
        self.create_text_input("stiffness_ds", input_frame, 'dsK0 [pN]', 3, 0, 2, 2)
        self.create_text_input("stiffness_ss", input_frame, 'ssK0 [pN]', 3, 2, 2, 2)
    
        ## handle part
        # persistance length for ds and ss
        self.create_text_input("dsLp", input_frame, 'dsLp [nm]', 4, 0, 2, 2)
        self.create_text_input("ssLp", input_frame, 'ssLp [nm]', 4, 2, 2, 2)
        # contour length for ds and ss
        self.create_text_input("dsLc", input_frame, 'dsLc [nm]', 5, 0, 2, 2)
        self.create_text_input("ssLc", input_frame, 'ssLc [nm]', 5, 2, 2, 2)


        # start and end positions of step
        label_strF = tk.Label(input_frame, text='Force [pN]', font='Helvetica 8 bold')
        label_strF.grid(row=6, column=1, sticky=tk.E + tk.W, pady=(25, 0))
        label_strD = tk.Label(input_frame, text='Distance [nm]', font='Helvetica 8 bold')
        label_strD.grid(row=6, column=3, sticky=tk.E + tk.W, pady=(25, 0))
        self.create_text_input("startF", input_frame, 'Start force [pN]', 7, 0, 2, 2)
        self.create_text_input("startD", input_frame, 'Start distance [nm]', 7, 2, 2, 2)
        self.create_text_input("endF", input_frame, 'End force [pN]', 8, 0, 2, 2)
        self.create_text_input("endD", input_frame, 'End distance [nm]', 8, 2, 2, 2)
        
        # create button widgets that use the defined functions
        self.create_button("Mark start", input_frame, self.start_click, 9, 0, 2, 0, 1, 10, 1)
        self.create_button("Mark end", input_frame, self.end_click, 10, 0, 2, 0, 1, 10, 1)
        self.create_button("Save step", input_frame, self.save_step, 9, 1, 2, 0, 1, 10, 1)
        self.create_button("Delete step", input_frame, self.delete_step, 10, 1, 2, 0, 1, 10, 1)

        # display help text for the shortcuts
        text_info = tk.Text(input_frame, height=10, width=30, background='#f0f0f0', wrap='word', borderwidth=0, highlightthickness=0)
        text_info.insert(tk.END, "Shortcuts", ('bold',))
        text_info.insert(tk.END, "\nmark step start")
        text_info.insert(tk.END, " <s>", ('bold',))
        text_info.insert(tk.END, "\nmark step end")
        text_info.insert(tk.END, " <e>", ('bold',))
        text_info.insert(tk.END, "\nsave step")
        text_info.insert(tk.END, " <Ctrl-s>", ('bold',))
        text_info.insert(tk.END, "\nstart analysis")
        text_info.insert(tk.END, " <Ctrl-f>", ('bold',))
        text_info.insert(tk.END, "\ndelete results line")
        text_info.insert(tk.END, " <mark+del>", ('bold',))
        text_info.insert(tk.END, "\nnext curve")
        text_info.insert(tk.END, " <Right arrow>", ('bold',))
        text_info.insert(tk.END, "\nprevious curve")
        text_info.insert(tk.END, " <Left arrow>", ('bold',))
        text_info.tag_configure('bold', font=('Arial', 10, 'bold'))
        text_info.grid(row=9, column=2, padx=4, pady=(20,20), columnspan=2, rowspan=2, sticky='nsew')
        text_info.configure(state='disabled')  # Make the text widget read-only

        # create Treeview for the steps to analyze
        cols_steps = ('Step number', 'F start', 'F end', 'Step start', 'Step end')
        tree_steps = ttk.Treeview(input_frame, columns=cols_steps, show='headings', height=5)

        # set column headings
        count = 0
        for col in cols_steps:
            tree_steps.heading(col, text=col)
            tree_steps.column(col, minwidth=25, width=65)

        tree_steps.grid(row=11, column=0, columnspan=2, pady=5)

        # step_number = 1
        # button_delete_step = tk.Button(TOMATO_parameter_frame, text='Delete step', command=delete_step, bg='lightsteelblue2', font=('Arial', 10, 'bold'))
        # button_delete_step.grid(row=8, column=0, pady=2)

        # button_start_analysis = tk.Button(TOMATO_parameter_frame, text='Analyze curve', command=analyze_steps, bg='#df4c4c', font=('Arial', 10, 'bold'))
        # button_start_analysis.grid(row=7, column=1, rowspan=2, pady=2)

        ### output frame
        self.create_text_input("entry_filename", TOMATO_results_frame, 'Filename', 0, 0, 4, 4, width=75)
 
        ## show the fitting parameters in a table
        # create Treeview for results table
        cols = (
            'Filename',
            'step number',
            'Force step start [pN]',
            'Force step end [pN]',
            'mean force [pN]',
            'extension step start [nm]',
            'extension step end [nm]',
            'Step length [nm]',
            'ds contour length',
            'ds persistance Length',
            'ds stiffness (K0) [pN]',
            'ss contour Length',
            'ss persistance Length',
            'ss stiffness (K0) [pN]',
            'Force offset',
            'Distance offset',
            'Work [pN*nm]',
            'Work [kT]'
        )

        tree_results = ttk.Treeview(TOMATO_results_frame, columns=cols, show='headings', height=5)
        # set column headings
        for col in cols:
            tree_results.heading(col, text=col)
            tree_results.column(col, minwidth=25, width=65)
        tree_results.grid(row=1, column=0, padx=5, pady=5, columnspan=4)


    def start_analysis(self):
        self.check_settings()


    def parameters(self, default_values, default_fit, default_constantF):
        if not default_values == 0:
            for key in self.text_input:
                if key in default_values:
                    self.text_input[key].delete(0, "end")
                    self.text_input[key].insert("end", str(default_values[key]))

        for key in self.text_input:
            if key in default_fit:
                self.text_input[key].delete(0, "end")
                self.text_input[key].insert("end", str(default_fit[key]))
            if key in default_constantF:
                self.text_input[key].delete(0, "end")
                self.text_input[key].insert("end", str(default_constantF[key]))


    def check_settings(self):
        wrong = 0
        for key in self.text_input:
            exclude_completly = ['startF', 'startD', 'endF', 'endD', 'entry_filename']
            exclude_non_numeric = ['mean_gauss', 'STD_gauss', 'amplitude_gauss']

            if not key in exclude_completly and not key in exclude_non_numeric: 
                parameter = self.text_input[key].get()

                try:
                    float(parameter)  # Attempt to convert parameter to float
                except ValueError:
                    self.print_output(f'Parameter {key} is not a number!')
                    wrong += 1
    
        if wrong == 0:
            return self.text_input, self.checkbox_vars
        else:
            return


    def create_check_box(self, key, frame, text, var_value, check_command, row, column, pad):
        check_var = tk.IntVar(value=var_value)
        check_box = tk.Checkbutton(frame, text=text, variable=check_var, command=check_command)
        check_box.grid(row=row, column=column, sticky='W', pady=pad, padx=pad)
        self.checkboxes[key] = check_box
        self.checkbox_vars[key] = check_var


    def create_text_input(self, key, frame, text, row, column, padx, pady, width=15):
        if not key in self.text_input.keys():
            text_var = tk.StringVar()
        else:
            text_var = self.text_input_vars[key]

        text_entry = tk.Entry(frame, textvariable=text_var, width=width)
        text_label = tk.Label(frame, text=text)
        text_label.grid(row=row, column=column, sticky='W', padx=padx, pady=pady)
        text_entry.grid(row=row, column=column+1, padx=padx, pady=pady)
        self.text_input[key] = text_entry
        self.text_input_vars[key] = text_var


    def create_button(self, text, frame, command, row, column, padx, pady, heigth, width, columnspan=1, cursor="hand2"):
        button = tk.Button(
            frame,
            text=text,
            command=command,
            bg='#df4c4c',
            activebackground='#eaa90d',
            font='Helvetica 10 bold',
            height=heigth,
            width=width,
            cursor=cursor
        )
        button.grid(row=row, column=column, padx=padx, pady=pady, columnspan=columnspan)


    def select_box(self, *selected_keys):
        for key1 in selected_keys:
            if self.checkbox_vars[key1].get() == 1:
                for key2 in selected_keys:
                    if not key1 == key2:
                        self.checkbox_vars[key2].set(value=0)

        boxes = [self.checkbox_vars[key].get() for key in selected_keys]
        if all(boxes) == 0:
            self.checkbox_vars[selected_keys[0]].set(value=1)


    def get_single_file(self):
        # Specify the allowed file types
        types = [
            ("Text or h5 files", "*.txt *.csv *.tsv *.h5"),
            ("All files", "*.*")
        ]

        import_file_path = tk.filedialog.askopenfilename(filetypes=types)
        file_extension = os.path.splitext(import_file_path)[1]

        if file_extension == '.csv' or file_extension == '.txt' or file_extension == '.tsv':
            if not self.checkbox_vars['CSV'].get() == 1:
                self.checkbox_vars['CSV'].set(value=1)
                self.select_box("CSV", "HF", "LF")
                # write message to output window
                self.print_output("CSV file selected for display, loaded CSV parameters.")
            elif self.checkbox_vars['CSV'].get() == 1:
                self.print_output("CSV file selected for display.")

                # parameters(default_values_CSV, default_values_FIT, default_values_constantF)

        elif file_extension == '.h5':
            if self.checkbox_vars['HF'].get() == 1:         
                # write message to output window
                self.print_output("H5 file selected (HF) for display.")
            elif self.checkbox_vars['LF'].get() == 1:
                # write message to output window
                self.print_output("H5 file selected (LF) for display.")

        text_params, checkbox_params = self.check_settings()

        checkbox_params['prepro'] = 0
        FD_raw, FD_raw_um, Frequency_value, filename = read_in_data(0, [import_file_path], text_params, checkbox_params)
        checkbox_params['prepro'] = 1
        FD, FD_um, Frequency_value, filename = read_in_data(0, [import_file_path], text_params, checkbox_params)
        # display_RAW_FD(FD[:, 0], FD[:, 1], FD_raw[:, 0], FD_raw[:, 1], filename)

    def print_output(self, text):
        current_time = time.localtime()
        current_time = time.strftime("%H:%M", current_time)
        self.output_window.config(state="normal")
        self.output_window.insert("end", f"{current_time} -- {text}\n")
        self.output_window.config(state="disabled")
        self.output_window.see("end")

    def refresh(self):
        self.root.update()
        self.root.update_idletasks()

    def show_h5(self):
        import_file_path = tk.filedialog.askopenfilename()
        h5_structure = lk.File(import_file_path)
        h5_structure_window = tk.Toplevel(self.root)
        h5_structure_window.title("H5 structure")

        text = tk.Text(h5_structure_window, height=50, width=200)
        scroll_bar = tk.Scrollbar(h5_structure_window, command=text.yview)
        scroll_bar.pack(side=tk.RIGHT, fill=tk.Y)
        text['yscrollcommand'] = scroll_bar.set
        text.pack(side=tk.LEFT, fill=tk.Y)
        text.insert("end", h5_structure)

    def load_parameters(self):
        pass

    def show_constantF(self):
        pass

    def start_constantF(self):
        pass

    def readme(self):
        pass

    def tab_bind(self, event):
        selected_tab = event.widget.nametowidget(event.widget.select())
        if selected_tab == self.tab_TOMATO:
            self.print_output("TOMATO is unresponsive during fitting, please be patient\n mark step start <s>\n mark step end <e>\n save marked step <Ctrl+s>\n start analysis <Ctrl+f>\n delete results line <mark+del>\n next curve <Right arrow>\n previous curve <Left arrow>")

    def show_augment(self):
        pass

    def start_TOMATO(self):
        pass

    def export_TOMATO(self):
        pass

    def reset_TOMATO(self):
        pass

    def start_click(self):
        pass

    def end_click(self):
        pass

    def save_step(self):
        pass
    
    def delete_step(self):
        pass

    def on_closing(self):
        # makes sure all python processes/loops are cancelled before exiting
        if tk.messagebox.askokcancel("Quit", "Do you really want to quit?"):
            self.root.destroy()


""" start the tkinter application """
if __name__ == '__main__':
    mp.freeze_support()
    app = FriedPotatoGUI()
    app.root.mainloop()





# input_settings = {
#     'downsample_value': int(downsample_value2.get()),
#     'filter_degree': int(Filter_degree2.get()),
#     'filter_cut_off': float(Filter_cut_off2.get()),
#     'F_min': float(Force_Min2.get()),
#     'step_d': int(step_d_value.get()),
#     'z-score_f': float(Z_score_force2.get()),
#     'z-score_d': float(Z_score_distance2.get()),
#     'window_size': int(window_size_value.get()),
#     'data_frequency': float(Frequency_value.get()),
#     'STD_diff': float(STD_difference_value.get()),
#     'augment_factor': augment_factor_value.get()
# }

# input_format = {
#     'HF': check_box_HF.get(),
#     'LF': check_box_LF.get(),
#     'CSV': check_box_CSV.get(),
#     'Augment': check_box_augment.get(),
#     'Trap': check_box_Trap1.get(),
#     'length_measure': check_box_um.get(),
#     'MultiH5': check_box_multiH5.get(),
#     'preprocess': check_box_preprocess.get(),
#     'reverse_fitting':check_box_reverse_fitting.get()
# }

# export_data = {
#     'export_SMOOTH': check_box_smooth_data.get(),
#     'export_PLOT': check_box_plot.get(),
#     'export_STEPS': check_box_steps.get(),
#     'export_TOTAL': check_box_total_results.get(),
#     'export_FIT': check_box_fitting.get()
# }

# input_fitting = {
#     'WLC+WLC': int(check_box_WLC.get()),
#     'WLC+FJC': int(check_box_FJC.get()),
#     'lp_ds': float(dsLp.get()),
#     'lp_ds_up': float(dsLp_up.get()),
#     'lp_ds_low': float(dsLp_low.get()),
#     'lc_ds': float(dsLc.get()),
#     'lp_ss': float(ssLp.get()),
#     'lc_ss': float(ssLc.get()),
#     'lc_ss_up': float(ssLc_up.get()),
#     'ds_stiff': float(stiff_ds.get()),
#     'ds_stiff_up': float(stiff_ds_up.get()),
#     'ds_stiff_low': float(stiff_ds_low.get()),
#     'ss_stiff': float(stiff_ss.get()),
#     'ss_stiff_up': float(stiff_ss_up.get()),
#     'ss_stiff_low': float(stiff_ss_low.get()),
#     'offset_f': float(f_off.get()),
#     'offset_f_up': float(f_off_up.get()),
#     'offset_f_low': float(f_off_low.get()),
#     'offset_d': float(d_off.get()),
#     'offset_d_up': float(d_off_up.get()),
#     'offset_d_low': float(d_off_low.get())
# }

# input_constantF = {
#     'x min': int(x_min.get()),
#     'x max': int(x_max.get()),
#     'y min': int(y_min.get()),
#     'y max': int(y_max.get()),
#     'Number gauss': int(number_gauss.get()),
#     'Mean': mean_gauss.get(),
#     'STD': STD_gauss.get(),
#     'Amplitude': amplitude_gauss.get()
# }

# return input_settings, input_format, export_data, input_fitting, input_constantF