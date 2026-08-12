"""FRIED_POTATO Results Browser - Interactive Data Visualization (Final Fixed v6)"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
import os
import glob
import lumicks.pylake as lk

class FRIED_POTATO_Browser:
    def __init__(self, root):
        self.root = root
        self.root.title("FRIED_POTATO Results Browser (v6)")
        self.root.geometry("1400x800")
        
        self.current_curve_idx = 0
        self.smooth_files = []
        self.total_results = None
        self.curve_data = {} # Stores {idx: {'force': ..., 'distance_nm': ..., 'filename': ...}}
        self.fitting_params = {} # Stores {idx: [list of ALL rows for this filename]}
        
        self._create_gui()
        
    def _create_gui(self):
        control_frame = tk.Frame(self.root)
        control_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.btn_open_folder = tk.Button(control_frame, text="Open Analysis Folder", command=self.open_analysis_folder,
            bg='#df4c4c', font='Helvetica 10 bold', height=2, width=20)
        self.btn_open_folder.pack(side=tk.LEFT, padx=5)
        
        self.status_label = tk.Label(control_frame, text="No folder loaded", fg='blue')
        self.status_label.pack(side=tk.LEFT, padx=20)
        
        nav_frame = tk.Frame(control_frame)
        nav_frame.pack(side=tk.RIGHT, padx=5)
        
        self.btn_prev = tk.Button(nav_frame, text="< Previous", command=self.prev_curve)
        self.btn_prev.pack(side=tk.LEFT, padx=2)
        
        self.lbl_curve_info = tk.Label(nav_frame, text="Curve 0 of 0", font='Helvetica 10 bold')
        self.lbl_curve_info.pack(side=tk.LEFT, padx=10)
        
        self.btn_next = tk.Button(nav_frame, text="Next >", command=self.next_curve)
        self.btn_next.pack(side=tk.LEFT, padx=2)
        
        main_container = tk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        left_panel = tk.Frame(main_container)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 5))
        
        params_tree_label = tk.Label(left_panel, text="Fitting Parameters", font='Helvetica 10 bold')
        params_tree_label.pack(anchor='w')
        
        columns = ('Parameter', 'Value', 'StdErr')
        self.params_tree = ttk.Treeview(left_panel, columns=columns, show='headings', height=15)
        for col in columns:
            self.params_tree.heading(col, text=col)
            self.params_tree.column(col, width=120)
        scrollbar = ttk.Scrollbar(left_panel, orient=tk.VERTICAL, command=self.params_tree.yview)
        self.params_tree.configure(yscrollcommand=scrollbar.set)
        self.params_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        right_panel = tk.Frame(main_container)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right_panel)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        toolbar = NavigationToolbar2Tk(self.canvas, right_panel)
        toolbar.update()
        
        bottom_frame = tk.Frame(self.root)
        bottom_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.btn_export_plot = tk.Button(bottom_frame, text="Export Current Plot", command=self.export_current_plot,
            bg='lightsteelblue2', font='Helvetica 9')
        self.btn_export_plot.pack(side=tk.LEFT, padx=5)
        
        self.btn_load_fits = tk.Button(bottom_frame, text="Reload Fit Parameters", command=self.load_all_curves,
            bg='palegreen2', font='Helvetica 9')
        self.btn_load_fits.pack(side=tk.LEFT, padx=5)
        
        self.root.bind('<Left>', lambda e: self.prev_curve())
        self.root.bind('<Right>', lambda e: self.next_curve())
        
    def open_analysis_folder(self):
        folder_path = filedialog.askdirectory(title="Select FRIED_POTATO Analysis Folder")
        if not folder_path: return
            
        try:
            total_results_files = glob.glob(os.path.join(folder_path, 'total_results_*.csv'))
            if not total_results_files:
                messagebox.showerror("Error", "No total_results file found!")
                return
            total_results_file = max(total_results_files, key=os.path.getctime)
            self.total_results = pd.read_csv(total_results_file)
            
            self.smooth_files = sorted(glob.glob(os.path.join(folder_path, '*_smooth_*.csv')))
            if not self.smooth_files:
                messagebox.showwarning("Warning", "No smooth files found!")
                
            self.load_all_curves()
            self.status_label.config(text=f"Loaded {len(self.smooth_files)} curves")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load analysis folder:\n{str(e)}")
            
    def load_all_curves(self):
        """
        Loads smooth files and groups ALL matching rows from total_results by filename.
        Ensures distance is converted from um -> nm.
        """
        self.curve_data.clear()
        self.fitting_params.clear()
        
        if self.total_results is None:
            return

        # Pre-group results by filename to avoid repeated filtering
        grouped_results = {}
        for _, row in self.total_results.iterrows():
            fname = str(row['filename'])
            if fname not in grouped_results:
                grouped_results[fname] = []
            grouped_results[fname].append(row.to_dict())

        for idx, smooth_file in enumerate(self.smooth_files):
            try:
                # Extract base_name by splitting on '_smooth_' and taking the part before it
                # e.g. "20250705-163335 FD Curve HD1_B01_P000_R03_T01_3r_smooth_20260616-180302.csv"
                #  -> "20250705-163335 FD Curve HD1_B01_P000_R03_T01_3r"
                basename_full = os.path.basename(smooth_file)  # get filename, drop path
                base_name = basename_full.split('_smooth_')[0]
                
                # Load Smooth Data and Convert Units (um -> nm)
                df = pd.read_csv(smooth_file)
                if len(df.columns) >= 2:
                    force = df.iloc[:, 0].values
                    distance_um = df.iloc[:, 1].values
                    
                    # CRITICAL: Multiply by 1000 to convert micrometers to nanometers
                    distance_nm = distance_um * 1000.0
                    
                    self.curve_data[idx] = {
                        'force': force, 
                        'distance': distance_nm,
                        'filename': base_name
                    }
                    
                    # Retrieve ALL fitting parameters for this file
                    matched_rows = grouped_results.get(base_name)
                    
                    if matched_rows:
                        self.fitting_params[idx] = matched_rows
                    else:
                        print(f"Info: No fitting results found for {base_name}")
                        
            except Exception as e:
                print(f"Warning: Failed to load {smooth_file}: {type(e).__name__}: {e}")
                
        self.update_navigation()
        if len(self.curve_data) > 0:
            self.display_current_curve()
            
    def update_navigation(self):
        total = len(self.curve_data)
        current = self.current_curve_idx + 1 if total > 0 else 0
        self.lbl_curve_info.config(text=f"Curve {current} of {total}")
        self.btn_prev.config(state=tk.NORMAL if total > 1 else tk.DISABLED)
        self.btn_next.config(state=tk.NORMAL if total > 1 else tk.DISABLED)
        
    def prev_curve(self):
        if self.current_curve_idx > 0:
            self.current_curve_idx -= 1
            self.display_current_curve()
            
    def next_curve(self):
        if self.current_curve_idx < len(self.curve_data) - 1:
            self.current_curve_idx += 1
            self.display_current_curve()
            
    def display_current_curve(self):
        if self.current_curve_idx not in self.curve_data:
            return
        self.ax.clear()
        
        data = self.curve_data[self.current_curve_idx]
        force = data['force']
        distance = data['distance'] # Already in nm
        
        # Plot the smooth curve
        self.ax.plot(distance, force, 'k-', linewidth=1, alpha=0.6, label='Smoothed FD-Curve')
        self.ax.set_xlabel('Distance (nm)', fontsize=12)
        self.ax.set_ylabel('Force (pN)', fontsize=12)
        self.ax.set_title(data['filename'], fontsize=14)
        self.ax.grid(True, alpha=0.3)
        
        if self.current_curve_idx in self.fitting_params:
            fit_params_list = self.fitting_params[self.current_curve_idx]
            if len(fit_params_list) > 0:
                # Pass the WHOLE list of rows (including ds and all steps)
                self._overlay_fits(fit_params_list, distance, force)
                # Update table with the FIRST row (ds part)
                self._update_params_table(fit_params_list[0])
                
        self.ax.legend(loc='best')
        self.fig.tight_layout()
        self.canvas.draw()
        
    def _safe_get_float(self, d, key, default=np.nan):
        val = d.get(key)
        if val is None or val == '' or val == 'nan' or val == 'NaN':
            return default
        try:
            f_val = float(val)
            return f_val if not np.isnan(f_val) else default
        except (ValueError, TypeError):
            return default

    def _overlay_fits(self, fit_params_list, distance, force):
        """
        Renders fits using corrected logic:
        1. ds_part uses ewlc_odijk_distance directly (better simulation behavior)
        2. Unfolded regions use the mixed DNA_2 + RNA model
        3. d_offset shifts the distance array before simulation
        4. Proper labeling for each region
        """
        if len(fit_params_list) == 0:
            return

        x_min = min(distance) - 100
        x_max = max(distance) + 100
        model_distance = np.linspace(x_min, x_max, 500)
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'cyan', 'magenta']
        
        # Extract ds_params from the FIRST row (index 0) - these are shared for ALL steps
        ds_params = fit_params_list[0]
        ds_Lc = self._safe_get_float(ds_params, 'Lc_ds', 100.0)
        ds_Lp = self._safe_get_float(ds_params, 'Lp_ds', 12.0)
        ds_St = self._safe_get_float(ds_params, 'St_ds', 750.0)
        ds_f_off = self._safe_get_float(ds_params, 'f_offset_ds', 0.0)
        ds_d_off = self._safe_get_float(ds_params, 'd_offset_ds', 0.0)
        
        for i, params in enumerate(fit_params_list[:10]):
            color = colors[i % len(colors)]
            model_type = str(params.get('model_type', 'UNKNOWN')).upper()
            
            try:
                # Helper to get values from current row
                def get_val(csv_key, default):
                    return self._safe_get_float(params, csv_key, default)
                
                sim_force = np.zeros_like(model_distance) * np.nan
                
                # --- CASE 1: ds_handle (First Row / i == 0) ---
                if i == 0:
                    # Use ewlc_odijk_distance directly (Distance -> Force)
                    # This avoids inversion issues
                    ds_model = lk.ewlc_odijk_distance("ds_part") + lk.force_offset("ds_part")
                    
                    fit_inst = lk.FdFit(ds_model)
                    # Add dummy data for infrastructure
                    fit_inst.add_data("sim", np.linspace(min(force), max(force), 5), model_distance[:5])
                    
                    # Set ds parameters
                    fit_inst["ds_part/Lc"].value = ds_Lc
                    fit_inst["ds_part/Lp"].value = ds_Lp
                    fit_inst["ds_part/St"].value = ds_St
                    fit_inst["ds_part/f_offset"].value = ds_f_off
                    
                    # Simulate (optionally apply d_offset shift)
                    effective_distance = model_distance - ds_d_off
                    sim_force = ds_model(effective_distance, fit_inst.params)
                    
                    label_name = "ds handle"
                    
                # --- CASE 2: Unfolded Regions (Rows 1+) ---
                elif "WLC+WLC" in model_type or "FULLY_UNFOLDED" in model_type:
                    # Determine if FJC or WLC for RNA
                    rna_model_str = "efjc_distance" if "FJC" in model_type else "ewlc_odijk_distance"
                    
                    # Build combined model: DNA_2 (ds from row 0) + RNA (from current row)
                    dna_model = getattr(lk, "ewlc_odijk_distance")("DNA_2")
                    rna_func = getattr(lk, rna_model_str)
                    rna_model = rna_func("RNA")
                    
                    combined = dna_model + rna_model
                    
                    # Apply inversion and offsets as requested
                    final_model = combined.invert().subtract_independent_offset() + lk.force_offset("DNA")
                    
                    fit_inst = lk.FdFit(final_model)
                    fit_inst.add_data("sim", np.linspace(min(force), max(force), 5), model_distance[:5])
                    
                    # Inject ds params from Row 0 into DNA_2 component
                    fit_inst["DNA_2/Lc"].value = ds_Lc
                    fit_inst["DNA_2/Lp"].value = ds_Lp
                    fit_inst["DNA_2/St"].value = ds_St
                    try:
                        fit_inst["DNA/f_offset"].value = ds_f_off
                    except KeyError: pass
                    
                    # Inject RNA params from CURRENT Row
                    fit_inst["RNA/Lc"].value = get_val('Lc_ss', 0.0)
                    fit_inst["RNA/Lp"].value = get_val('Lp_ss', 1.0)
                    fit_inst["RNA/St"].value = get_val('St_ss', 1000.0)
                    
                    # Get d_offset from current row (if different from ds)
                    step_d_off = get_val('d_offset_ds', ds_d_off)
                    
                    # Simulate with distance offset shift
                    effective_distance = model_distance - step_d_off
                    
                    try:
                        sim_force = final_model(effective_distance, fit_inst.params)
                    except:
                        # If simulation fails, skip this fit
                        continue
                    
                    # Label based on index and model type
                    if "FULLY_UNFOLDED" in model_type:
                        label_name = "Fully Unfolded"
                    else:
                        label_name = f"Step {i}"
                
                elif "WLC+FJC" in model_type:
                    # Similar to WLC+WLC but with FJC for RNA
                    rna_func = getattr(lk, "efjc_distance")
                    
                    dna_model = getattr(lk, "ewlc_odijk_distance")("DNA_2")
                    rna_model = rna_func("RNA")
                    
                    combined = dna_model + rna_model
                    final_model = combined.invert().subtract_independent_offset() + lk.force_offset("DNA")
                    
                    fit_inst = lk.FdFit(final_model)
                    fit_inst.add_data("sim", np.linspace(min(force), max(force), 5), model_distance[:5])
                    
                    fit_inst["DNA_2/Lc"].value = ds_Lc
                    fit_inst["DNA_2/Lp"].value = ds_Lp
                    fit_inst["DNA_2/St"].value = ds_St
                    try:
                        fit_inst["DNA/f_offset"].value = ds_f_off
                    except KeyError: pass
                    
                    fit_inst["RNA/Lc"].value = get_val('Lc_ss', 0.0)
                    fit_inst["RNA/Lp"].value = get_val('Lp_ss', 1.0)
                    fit_inst["RNA/St"].value = get_val('St_ss', 1000.0)
                    
                    step_d_off = get_val('d_offset_ds', ds_d_off)
                    effective_distance = model_distance - step_d_off
                    
                    try:
                        sim_force = final_model(effective_distance, fit_inst.params)
                    except:
                        continue
                        
                    label_name = f"Step {i} (FJC)"
                
                else:
                    continue
                
                # Filter valid points
                valid_mask = np.isfinite(sim_force) & (model_distance >= x_min) & (model_distance <= x_max)
                if np.any(valid_mask):
                    self.ax.plot(model_distance[valid_mask], sim_force[valid_mask], 
                               linestyle='--', color=color, linewidth=2, alpha=0.8, label=label_name)
                               
            except Exception as e:
                print(f"Debug: Fit {i} failed ({model_type}): {e}")
                pass
            
    def _update_params_table(self, params_dict):
        for item in self.params_tree.get_children():
            self.params_tree.delete(item)
        key_params = [
            ('Model Type', params_dict.get('model_type', 'N/A')),
            ('Log Likelihood', params_dict.get('log_likelihood', 'N/A')),
            ('Lc_ds (nm)', params_dict.get('Lc_ds', 'N/A')),
            ('Lp_ds (nm)', params_dict.get('Lp_ds', 'N/A')),
            ('St_ds (pN)', params_dict.get('St_ds', 'N/A')),
            ('Lc_ss (nm)', params_dict.get('Lc_ss', 'N/A')),
            ('Lp_ss (nm)', params_dict.get('Lp_ss', 'N/A')),
            ('St_ss (pN)', params_dict.get('St_ss', 'N/A')),
            ('Work (kT)', params_dict.get('Work_(kB*T)', 'N/A')),
        ]
        for param, value in key_params:
            std_val = "" # Std error column currently empty
            self.params_tree.insert('', 'end', values=(param, str(value), str(std_val)))
            
    def export_current_plot(self):
        save_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png"), ("SVG files", "*.svg"), ("All files", "*.*")])
        if save_path:
            try:
                if save_path.endswith('.svg'):
                    self.fig.savefig(save_path, format='svg', dpi=600)
                else:
                    self.fig.savefig(save_path, dpi=600)
                messagebox.showinfo("Success", f"Plot saved to:\n{save_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save plot:\n{e}")

def main():
    root = tk.Tk()
    app = FRIED_POTATO_Browser(root)
    root.protocol("WM_DELETE_WINDOW", lambda: exit_app(app))
    root.mainloop()

def exit_app(app):
    app.root.destroy()
    import sys
    sys.exit(0)

if __name__ == "__main__":
    main()