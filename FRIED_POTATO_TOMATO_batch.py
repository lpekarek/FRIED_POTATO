"""Copyright 2026 Lukáš Pekárek & Stefan Buck & Simon Doll & Leo König"""

"""Headless TOMATO batch worker.

Runs 'Analyse all marked curves' in a separate process, mirroring the classic
FRIED_POTATO_ForceRamp.start_subprocess() pattern: no tkinter, progress sent
via output_q, all outputs written into the Manual_analysis_<timestamp> folder.
This is a faithful headless port of the GUI's analyze_steps() fitting loop.
"""

import os
import time
import traceback
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # headless: the subprocess must never open a window
from matplotlib.figure import Figure

from FRIED_POTATO_ForceRamp import read_in_data_TOMATO
from FRIED_POTATO_preprocessing import create_derivative
from FRIED_POTATO_find_steps import calc_integral
import FRIED_POTATO_fitting
from FRIED_POTATO_fitting import fitting_ds, fitting_ss, fitting_FU, fitting_FU_ss


TOMATO_TOTAL_HEADER = [
    'filename', 'orientation', 'Derivative of', 'step number',
    'F1', 'F2', 'Fc', 'step start', 'step end', 'step length',
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


def log_error(error_msg, analysis_folder, timestamp):
    """Self-contained error log (the GUI's log_error isn't importable)."""
    try:
        with open(os.path.join(analysis_folder, f'error_log_{timestamp}.txt'), 'a') as f:
            f.write(error_msg + '\n')
    except Exception:
        pass


def run_tomato_batch(job, output_q):
    """Entry point of the subprocess. job is a plain dict (picklable)."""
    analysis_folder = job['analysis_folder']
    timestamp = job['timestamp']
    Files = job['files']
    marks = job['marks']
    input_settings = job['input_settings']
    input_format = job['input_format']
    export_data = job['export_data']
    input_fitting = job['input_fitting']
    sort_steps = job.get('sort_steps', False)

    # cumulative total_results file with the batch-compatible header
    filename_total_results = os.path.join(
        analysis_folder, 'total_results_' + timestamp + '.csv')
    with open(filename_total_results, 'w') as f:
        f.write(','.join(TOMATO_TOTAL_HEADER))
        f.write('\n')

    entries = sorted(marks.items(), key=lambda kv: kv[1]['index'])
    n_total = len(entries)
    output_q.put(f'TOMATO batch started: {n_total} curve(s), '
                 f'output in {os.path.basename(analysis_folder)}')

    for k, (fname_hint, rec) in enumerate(entries):
        try:
            # --- load the curve exactly like change_FD()/open_folder() do ---
            FD, FD_um, Freq, filename_i, FD_ds = read_in_data_TOMATO(
                rec['index'], Files, input_settings, input_format)
            der_arr = create_derivative(
                input_settings, Freq, FD[:, 0], FD[:, 1], 0)

            orientation_i = 'forward'
            if len(FD) >= 2 and FD[0, 1] > FD[-1, 1]:
                FD = np.flipud(FD)
                orientation_i = 'reverse'

            output_q.put(f'TOMATO: fitting {filename_i} ({k + 1}/{n_total})')

            export_df = analyze_one_curve(
                filename_i, FD, der_arr, orientation_i, rec['rows'],
                input_settings, input_format, export_data, input_fitting,
                analysis_folder, timestamp, sort_steps)

            export_df.to_csv(filename_total_results, mode='a',
                             index=False, header=False)

        except Exception as e:
            msg = (f'ERROR while fitting {fname_hint}: {e}\n'
                   f'{traceback.format_exc()}')
            print(msg)                     # <-- terminal output
            output_q.put(f'TOMATO ERROR: {fname_hint} failed - {e}')
            log_error(msg, analysis_folder, timestamp)

    output_q.put('TOMATO batch fitting finished. '
                 f'Results in: {analysis_folder}')


def analyze_one_curve(filename_i, FD, der_arr, orientation_i, steps,
                      input_settings, input_format, export_data, input_fitting,
                      out_dir, timestamp, sort_steps):
    """Headless port of analyze_steps(): fits every region defined by the
    marked steps, computes Work, writes per-curve fit CSV / PNG (+SVG) /
    pickle into out_dir, and returns the batch-compatible export DataFrame."""

    if not steps:
        raise ValueError('No steps marked for this curve.')

    # ---- step table (replaces the tree_steps treeview) ----------------------
    columns = ('Step number', 'F start', 'F end', 'Step start', 'Step end')
    # rows may arrive as 5-tuples (Step number, F start, F end, Step start,
    # Step end) or 4-tuples (step number stripped by the sender) - normalise:
    # always take the LAST four values and derive a consecutive step number.
    data_rows = []
    for n, r in enumerate(steps, start=1):
        vals = list(r)[-4:]              # F start, F end, Step start, Step end
        data_rows.append((n,) + tuple(float(v) for v in vals))
    treeview_df = pd.DataFrame(data_rows, columns=columns)

    # optional sorting by 'Step end' (same as the GUI checkbox)
    if sort_steps:
        ends = pd.to_numeric(treeview_df['Step end'], errors='coerce')
        treeview_df = (treeview_df.assign(_end_num=ends)
                       .sort_values('_end_num')
                       .drop(columns='_end_num')
                       .reset_index(drop=True))
        treeview_df['Step number'] = range(1, len(treeview_df) + 1)

    # ---- figure (Figure object only - savefig works without any GUI) --------
    figure1 = Figure(figsize=(8, 6), dpi=110)
    diff_colors = ['b', 'r', 'c', 'g', 'y', 'm'] * 4
    subplot1 = figure1.add_subplot(111)
    subplot1.plot(FD[:, 1], FD[:, 0], color='gray')
    distance = np.arange(min(FD[:, 1]), max(FD[:, 1]) + 50, 2)
    fit_data = {"distance": distance}
    export_fit = []
    fit = []
    start_force_ss = []
    start_distance_ss = []
    integral_ss_fit_start = []
    integral_ss_fit_end = []
    rows_for_export = []

    def make_export_row(fit_dict, step_number, F1='', F2='', Fc='',
                        step_start='', step_end='', step_length=''):
        return {
            'filename': filename_i,
            'orientation': orientation_i,
            'Derivative of': '',
            'step number': step_number,
            'F1': F1, 'F2': F2, 'Fc': Fc,
            'step start': step_start, 'step end': step_end,
            'step length': step_length,
            'model_type': fit_dict.get('model_type', ''),
            'log_likelihood': fit_dict.get('log_likelihood', ''),
            'Lc_ds': fit_dict.get('Lc_ds', ''),
            'Lc_ds_stderr': fit_dict.get('Lc_ds_stderr', ''),
            'Lp_ds': fit_dict.get('Lp_ds', ''),
            'Lp_ds_stderr': fit_dict.get('Lp_ds_stderr', ''),
            'St_ds': fit_dict.get('St_ds', ''),
            'St_ds_stderr': fit_dict.get('St_ds_stderr', ''),
            'f_offset_ds': fit_dict.get('f_offset_ds', ''),
            'f_offset_ds_stderr': fit_dict.get('f_offset_ds_stderr', ''),
            'd_offset_ds': fit_dict.get('d_offset_ds', ''),
            'd_offset_ds_stderr': fit_dict.get('d_offset_ds_stderr', ''),
            'Lc_ss': fit_dict.get('Lc_ss', ''),
            'Lc_ss_stderr': fit_dict.get('Lc_ss_stderr', ''),
            'Lp_ss': fit_dict.get('Lp_ss', ''),
            'Lp_ss_stderr': fit_dict.get('Lp_ss_stderr', ''),
            'St_ss': fit_dict.get('St_ss', ''),
            'St_ss_stderr': fit_dict.get('St_ss_stderr', ''),
            'f_offset_ss': fit_dict.get('f_offset_ss', ''),
            'f_offset_ss_stderr': fit_dict.get('f_offset_ss_stderr', ''),
            'd_offset_ss': fit_dict.get('d_offset_ss', ''),
            'd_offset_ss_stderr': fit_dict.get('d_offset_ss_stderr', ''),
            'Work_(pN*nm)': '', 'Work_(kB*T)': '',
            'delta Lc': '', 'total Lc': '', 'total W': '',
            'total number of steps': '',
            'fit_status': fit_dict.get('fit_status', '')
        }

    reverse_fitting_enabled = input_format.get('reverse_fitting', 0) == 1

    for i in treeview_df.index:
        if treeview_df['Step number'][i] == 1:
            j = treeview_df['Step number'][i]

            if reverse_fitting_enabled:
                ds_fit_dict, TOMATO_area_ds, real_step_start, model_FU_obj = fitting_FU(
                    filename_i, input_settings, export_data, input_fitting,
                    float(treeview_df['Step start'][i]),
                    float(treeview_df['Step end'][len(treeview_df) - 1]),
                    max(FD[:, 1]), FD, der_arr, [], 1)

                fit_ss, f_fitting_region_ss, d_fitting_region_ss, \
                    first_region_dict, area_ss_fit_start, area_ss_fit_end = fitting_FU_ss(
                        filename_i, input_settings, export_data, input_fitting,
                        0, float(treeview_df['Step start'][i]),
                        FD, 1, 1, der_arr, [], 1)

                fit.append(fit_ss)
                start_force_ss.append(f_fitting_region_ss)
                start_distance_ss.append(d_fitting_region_ss)
                integral_ss_fit_start.append(area_ss_fit_start)
                integral_ss_fit_end.append(area_ss_fit_end)
                export_fit.append(first_region_dict)
                rows_for_export.append(make_export_row(first_region_dict, step_number=0))


                F_ss_model = FRIED_POTATO_fitting.model_ss(distance, fit_ss.params)
                fit_data[f"Fit_ss_{i}"] = F_ss_model
                subplot1.plot(d_fitting_region_ss, f_fitting_region_ss, color=diff_colors[i])
                subplot1.plot(distance, F_ss_model, marker=None, linewidth=1, linestyle='dashed', color="black")

                # --- NEW: plot the ds model curve from the FU fit, dashed grey, on top ---
                try:
                    F_ds_model = FRIED_POTATO_fitting.model_ds(
                        distance, FRIED_POTATO_fitting.fit_ds.params)
                    fit_data[f"Fit_ds_{i}"] = F_ds_model
                    subplot1.plot(distance, F_ds_model, marker=None, linewidth=1.5,
                                linestyle='dashed', color='grey', zorder=10)
                except Exception as e:
                    print(f"Warning: could not plot ds model curve: {e}")

                ds_fit_region_end = real_step_start

            else:
                ds_fit_dict, TOMATO_area_ds, real_step_start, fit_ds = fitting_ds(
                    filename_i, input_settings, export_data, input_fitting,
                    float(treeview_df['Step start'][i]), FD, der_arr, [], 1)

                F_ds_model = FRIED_POTATO_fitting.model_ds(distance, fit_ds.params)
                subplot1.plot(FD[:, 1][:real_step_start], FD[:, 0][:real_step_start],
                              color=diff_colors[i])
                subplot1.plot(distance, F_ds_model, linestyle='dashed',
                              linewidth=1, color='black')
                ds_fit_region_end = real_step_start

            if not reverse_fitting_enabled:
                rows_for_export.append(make_export_row(ds_fit_dict, step_number=0))
                export_fit.append(ds_fit_dict)

        elif treeview_df['Step number'][i] > 1:
            j = treeview_df['Step number'][i]

            fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, \
                area_ss_fit_start, area_ss_fit_end = fitting_ss(
                    filename_i, input_settings, export_data, input_fitting,
                    float(treeview_df['Step end'][i - 1]),
                    float(treeview_df['Step start'][i]),
                    FD, 1, 1, der_arr, [], 1)

            fit.append(fit_ss)
            start_force_ss.append(f_fitting_region_ss)
            start_distance_ss.append(d_fitting_region_ss)
            export_fit.append(ss_fit_dict)
            integral_ss_fit_start.append(area_ss_fit_start)
            integral_ss_fit_end.append(area_ss_fit_end)

            rows_for_export.append(make_export_row(
                ss_fit_dict, step_number=j - 1,
                F1=FD[real_step_start, 0],
                F2=f_fitting_region_ss[0],
                Fc=(f_fitting_region_ss[0] + FD[real_step_start, 0]) / 2,
                step_start=FD[real_step_start, 1],
                step_end=d_fitting_region_ss[0],
                step_length=d_fitting_region_ss[0] - FD[real_step_start, 1]))

            real_step_start = np.where(FD[:, 0] == f_fitting_region_ss[-1])[0][0]

            F_ss_model = FRIED_POTATO_fitting.model_ss(distance, fit_ss.params)
            fit_data[f"Fit_ss_{i}"] = F_ss_model
            subplot1.plot(d_fitting_region_ss[:], f_fitting_region_ss, color=diff_colors[i])
            subplot1.plot(distance, F_ss_model, linewidth=1, linestyle='dashed', color='black')

    # ---- last region of the curve -------------------------------------------
    fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, \
        area_ss_fit_start, area_ss_fit_end = fitting_ss(
            filename_i, input_settings, export_data, input_fitting,
            float(treeview_df['Step end'][len(treeview_df) - 1]),
            max(FD[:, 1]), FD, 1, 1, der_arr, [], 1)

    fit.append(fit_ss)
    start_force_ss.append(f_fitting_region_ss)
    start_distance_ss.append(d_fitting_region_ss)
    export_fit.append(ss_fit_dict)
    integral_ss_fit_start.append(area_ss_fit_start)
    integral_ss_fit_end.append(area_ss_fit_end)

    rows_for_export.append(make_export_row(
        ss_fit_dict, step_number=j,
        F1=FD[:, 0][real_step_start],
        F2=f_fitting_region_ss[0],
        Fc=(f_fitting_region_ss[0] + FD[:, 0][real_step_start]) / 2,
        step_start=FD[:, 1][real_step_start],
        step_end=d_fitting_region_ss[0],
        step_length=d_fitting_region_ss[0] - FD[:, 1][real_step_start]))

    F_ss_model = FRIED_POTATO_fitting.model_ss(distance, fit_ss.params)
    subplot1.plot(d_fitting_region_ss[:], f_fitting_region_ss, color=diff_colors[j + 1])
    subplot1.plot(distance, F_ss_model, linewidth=1, linestyle='dashed', color='black')

    # ---- work per step -------------------------------------------------------
    work_first_step, kT_1 = calc_integral(
        TOMATO_area_ds, integral_ss_fit_start[0],
        FD[ds_fit_region_end, 1], start_distance_ss[0][0],
        FD[ds_fit_region_end, 0], start_force_ss[0][0])

    for row in rows_for_export:
        if row['step number'] == 1:
            row['Work_(pN*nm)'] = work_first_step
            row['Work_(kB*T)'] = kT_1

    if j > 1:
        for n in range(1, j):
            work_step_n, kT_n = calc_integral(
                integral_ss_fit_end[n - 1], integral_ss_fit_start[n],
                start_distance_ss[n - 1][-1], start_distance_ss[n][0],
                start_force_ss[n - 1][-1], start_force_ss[n][0])
            for row in rows_for_export:
                if row['step number'] == n + 1:
                    row['Work_(pN*nm)'] = work_step_n
                    row['Work_(kB*T)'] = kT_n

    # ---- batch-style derived fields ------------------------------------------
    export_df = pd.DataFrame(rows_for_export, columns=TOMATO_TOTAL_HEADER)
    for col in ['Lc_ss', 'Work_(pN*nm)', 'Work_(kB*T)']:
        export_df[col] = pd.to_numeric(export_df[col], errors='coerce')

    export_df['delta Lc'] = export_df['Lc_ss'].diff().fillna("#N/A")
    export_df['total Lc'] = "#N/A"
    export_df['total W'] = "#N/A"
    export_df['total number of steps'] = "#N/A"
    try:
        export_df.loc[export_df.index[-1], 'total Lc'] = export_df['Lc_ss'].iloc[-1]
    except (IndexError, KeyError):
        pass
    try:
        total_W = export_df['Work_(kB*T)'].sum()
        export_df.loc[export_df.index[-1], 'total W'] = \
            "#N/A" if pd.isna(total_W) else total_W
    except (IndexError, KeyError):
        pass
    try:
        export_df.loc[export_df.index[-1], 'total number of steps'] = \
            export_df['step number'].iloc[-1]
    except (IndexError, KeyError):
        pass
    export_df = export_df.replace({np.nan: ''})

    # ---- per-curve exports (fit CSV, PNG/SVG, pickle) ------------------------
    fit_df = pd.DataFrame(fit_data)
    fit_df.to_csv(os.path.join(out_dir, f"{filename_i}_fit_data_{timestamp}.csv"),
                  index=False)

    subplot1.set_ylim([min(FD[:, 0]), max(FD[:, 0])])
    subplot1.set_xlim([min(FD[:, 1]) - 10, max(FD[:, 1]) + 10])
    subplot1.tick_params('both', direction='in')
    subplot1.set_xlabel("Distance [nm]")
    subplot1.set_ylabel("Force [pN]")

    figure1.savefig(os.path.join(out_dir, f"{filename_i}_fit_{timestamp}.png"), dpi=300)
    if export_data.get('export_svg', 0) == 1:
        figure1.savefig(os.path.join(out_dir, f"{filename_i}_fit_{timestamp}.svg"), dpi=300)

    with open(os.path.join(out_dir, f"{filename_i}_fit_objects_{timestamp}.pickle"),
              'wb') as f:
        pickle.dump({
            'fit_objects': fit,
            'fit_dicts': export_fit,
            'filename': filename_i,
            'steps': treeview_df
        }, f)

    return export_df