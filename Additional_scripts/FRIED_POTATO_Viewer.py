import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import glob
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import lumicks.pylake as lk

class FriedPotatoViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("FRIED POTATO - Interactive Viewer & Model Explorer")
        self.root.geometry("1400x800")

        # State variables
        self.analysis_folder = None
        self.total_df = None
        self.smooth_files = {}
        self.pkl_files = {}
        self.current_filename = None
        self.current_data = None
        self.current_models = None
        
        self._create_widgets()
        
    def _create_widgets(self):
        top_frame = tk.Frame(self.root, pady=10)
        top_frame.pack(side=tk.TOP, fill=tk.X)
        
        btn_browse = tk.Button(top_frame, text="Select Analysis Folder", command=self.load_folder, 
                               bg="#df4c4c", fg="white", font=('Arial', 10, 'bold'))
        btn_browse.pack(side=tk.LEFT, padx=20)
        
        self.lbl_status = tk.Label(top_frame, text="No folder selected", font=('Arial', 10))
        self.lbl_status.pack(side=tk.LEFT, padx=20)

        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left Pane: Plot
        left_frame = tk.Frame(main_paned, width=900)
        main_paned.add(left_frame, weight=3)
        
        self.plot_fig = Figure(figsize=(10, 6), dpi=100)
        self.plot_ax = self.plot_fig.add_subplot(111)
        self.plot_canvas = FigureCanvasTkAgg(self.plot_fig, master=left_frame)
        self.plot_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        toolbar_frame = tk.Frame(left_frame)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.plot_canvas, toolbar_frame)
        self.toolbar.update()
        
        nav_frame = tk.Frame(left_frame, pady=5)
        nav_frame.pack(side=tk.BOTTOM, fill=tk.X)
        
        btn_prev = tk.Button(nav_frame, text="<< Previous Curve", command=self.prev_curve)
        btn_prev.pack(side=tk.LEFT, padx=5)
        
        lbl_curr = tk.Label(nav_frame, text="No curve loaded", font=('Arial', 9))
        self.lbl_curve_name = lbl_curr
        lbl_curr.pack(side=tk.LEFT, padx=10)
        
        btn_next = tk.Button(nav_frame, text="Next Curve >>", command=self.next_curve)
        btn_next.pack(side=tk.RIGHT, padx=5)

        # Right Pane: Table
        right_frame = tk.Frame(main_paned, width=400)
        main_paned.add(right_frame, weight=2)
        
        tree_frame = tk.Frame(right_frame)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        cols = ('Step #', 'F1 (pN)', 'F2 (pN)', 'Start (nm)', 'End (nm)', 
                'Lc_ds', 'Lp_ds', 'St_ds', 'f_off', 'd_off', 
                'Lc_ss', 'Lp_ss', 'St_ss')
        
        self.tree = ttk.Treeview(tree_frame, columns=cols, show='headings')
        for col in cols:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=70, anchor=tk.CENTER)
            
        v_scroll = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        h_scroll = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=v_scroll.set, xscrollcommand=h_scroll.set)
        
        v_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        h_scroll.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.tree.bind('<<TreeviewSelect>>', self.on_step_select)

        info_frame = tk.Frame(right_frame, pady=5)
        info_frame.pack(fill=tk.X)
        
        self.info_text = tk.Text(info_frame, height=8, wrap=tk.WORD, font=('Courier', 9))
        self.info_text.pack(fill=tk.X)
        
    def load_folder(self):
        folder = filedialog.askdirectory(title="Select FRIED POTATO Analysis Folder")
        if not folder:
            return
            
        self.analysis_folder = folder
        self.lbl_status.config(text=f"Folder: {os.path.basename(folder)}")
        
        total_file = glob.glob(os.path.join(folder, "total_results_*.csv"))
        if not total_file:
            messagebox.showerror("Error", "No 'total_results_*.csv' found!")
            return
            
        try:
            self.total_df = pd.read_csv(total_file[0])
            print(f"Loaded {len(self.total_df)} rows.")
        except Exception as e:
            messagebox.showerror("Error", f"Could not read CSV: {e}")
            return
        
        self.smooth_files = {}
        self.pkl_files = {}
        
        unique_files = self.total_df['filename'].unique()
        
        for fname in unique_files:
            smooth_candidates = glob.glob(os.path.join(folder, f"{fname}_smooth_*.csv"))
            pkl_candidates = glob.glob(os.path.join(folder, f"{fname}_pylake_models_*.pkl"))
            
            if not smooth_candidates:
                all_smooths = glob.glob(os.path.join(folder, "*_smooth_*.csv"))
                for s in all_smooths:
                    if fname.replace('_fw_curve', '').replace('_rv_curve', '') in os.path.basename(s):
                        smooth_candidates = [s]
                        break
                        
            if not pkl_candidates:
                all_pkls = glob.glob(os.path.join(folder, "*_pylake_models_*.pkl"))
                for p in all_pkls:
                    if fname.replace('_fw_curve', '').replace('_rv_curve', '') in os.path.basename(p):
                        pkl_candidates = [p]
                        break
            
            if smooth_candidates:
                self.smooth_files[fname] = smooth_candidates[0]
            if pkl_candidates:
                self.pkl_files[fname] = pkl_candidates[0]
        
        self.curve_list = list(unique_files)
        self.current_idx = 0
        
        if self.curve_list:
            self.load_curve(0)
        else:
            messagebox.showwarning("Empty", "No curves found.")

    def load_curve(self, idx):
        if not self.curve_list or idx < 0 or idx >= len(self.curve_list):
            return
            
        self.current_idx = idx
        self.current_filename = self.curve_list[idx]
        self.lbl_curve_name.config(text=f"Curve: {self.current_filename}")
        
        # Load Smooth Data
        smooth_path = self.smooth_files.get(self.current_filename)
        self.current_data = None
        if smooth_path and os.path.exists(smooth_path):
            try:
                # Assuming CSV has NO header based on generation logic in FRIED_POTATO_ForceRamp.py
                self.current_data = pd.read_csv(smooth_path, header=None)
                print(f"Loaded smooth data: shape {self.current_data.shape}")
            except Exception as e:
                print(f"Error loading smooth: {e}")
                self.current_data = None

        # Load Models
        pkl_path = self.pkl_files.get(self.current_filename)
        self.current_models = None
        if pkl_path and os.path.exists(pkl_path):
            try:
                with open(pkl_path, 'rb') as f:
                    self.current_models = pickle.load(f)
                print(f"Loaded models: {list(self.current_models.keys())}")
            except Exception as e:
                print(f"Error loading PKL: {e}")
                self.current_models = None

        # Populate Table
        curve_rows = self.total_df[self.total_df['filename'] == self.current_filename]
        self.tree.delete(*self.tree.get_children())
        
        for _, row in curve_rows.iterrows():
            def safe_float(val, default=-999):
                try:
                    v = float(val)
                    return v if not np.isnan(v) else default
                except: return default

            values = [
                int(row.get('step number', 0)),
                safe_float(row.get('F1')), safe_float(row.get('F2')),
                safe_float(row.get('step start')), safe_float(row.get('step end')),
                safe_float(row.get('Lc_ds')), safe_float(row.get('Lp_ds')), safe_float(row.get('St_ds')),
                safe_float(row.get('f_offset_ds')), safe_float(row.get('d_offset_ds')),
                safe_float(row.get('Lc_ss')), safe_float(row.get('Lp_ss')), safe_float(row.get('St_ss'))
            ]
            
            disp_values = ["-" if v == -999 else f"{v:.2f}" for v in values]
            self.tree.insert("", "end", values=tuple(disp_values))
            
        if self.tree.get_children():
            self.tree.selection_set(self.tree.get_children()[0])
            self.on_step_select(None)

    def on_step_select(self, event):
        if self.current_data is None or (hasattr(self.current_data, 'empty') and self.current_data.empty):
            self.info_text.insert(tk.END, "No data available.\n")
            return
        if self.current_models is None:
            self.info_text.insert(tk.END, "No models loaded.\n")
            return
        
        selection = self.tree.selection()
        if not selection:
            return
            
        curve_rows = self.total_df[self.total_df['filename'] == self.current_filename].reset_index(drop=True)
        selected_row_idx = self.tree.index(selection[0])
        if selected_row_idx >= len(curve_rows):
            return

        row = curve_rows.iloc[selected_row_idx]
        step_num = int(row['step number'])
        
        self.info_text.delete("1.0", tk.END)
        self.info_text.insert(tk.END, f"Step {step_num}: {row['model_type']}\n")
        
        f1, f2 = row.get('F1', np.nan), row.get('F2', np.nan)
        start_dist, end_dist = row.get('step start', np.nan), row.get('step end', np.nan)
        
        self.info_text.insert(tk.END, f"F1: {f1:.2f}, F2: {f2:.2f} | Range: {start_dist:.2f}-{end_dist:.2f} nm\n")

        self.plot_ax.clear()
        self.plot_ax.set_title(f"{self.current_filename} - Step {step_num}", fontsize=14)
        self.plot_ax.set_xlabel("Distance (nm)")
        self.plot_ax.set_ylabel("Force (pN)")

        D_vals = self.current_data.iloc[:, 1] if self.current_data.shape[1] > 1 else None
        F_vals = self.current_data.iloc[:, 0] if self.current_data.shape[1] > 0 else None
        
        if D_vals is not None and F_vals is not None:
            mask = ~np.isnan(F_vals) & ~np.isnan(D_vals)
            self.plot_ax.scatter(D_vals[mask], F_vals[mask], s=5, alpha=0.3, color='gray', label='Data')

        # Plot ALL fits from the pickle file
        colors = ['b', 'r', 'c', 'g', 'm', 'y', 'k']
        fit_count = 0
        
        for key, fit_obj in self.current_models.items():
            if key.startswith('_'): continue # Skip internal keys
            
            try:
                # Determine min/max range from data
                min_d, max_d = (D_vals.min(), D_vals.max()) if D_vals is not None else (0, 1000)
                dist_plot = np.linspace(min_d, max_d, 500)
                
                force_plot = None
                
                # Try to evaluate
                try:
                    # Method 1: Direct call (inverted model returns Force for Distance)
                    force_plot = fit_obj(dist_plot)
                except Exception as e1:
                    try:
                        # Method 2: Access via .models
                        if hasattr(fit_obj, 'models'):
                            m = list(fit_obj.models.values())[0] if isinstance(fit_obj.models, dict) else fit_obj.models[0]
                            force_plot = m(dist_plot, fit_obj.params)
                        else:
                            raise
                    except Exception as e2:
                        try:
                            # Method 3: .model attribute
                            force_plot = fit_obj.model(dist_plot, fit_obj.params)
                        except:
                            pass
                
                if force_plot is not None:
                    color = colors[fit_count % len(colors)]
                    self.plot_ax.plot(dist_plot, force_plot, linestyle='--', color=color, linewidth=1.5, label=f"{key}")
                    
                    # Highlight region if this key corresponds to the selected step
                    # We assume ss_fit_N corresponds to step N+1
                    target_key = f"ss_fit_{step_num-1}" if step_num > 0 else "ds_initial_fit"
                    if key == target_key and not np.isnan(start_dist) and not np.isnan(end_dist):
                        self.plot_ax.axvspan(float(start_dist), float(end_dist), alpha=0.15, color=color, label='Selected Region')
                        self.plot_ax.legend(loc='upper left', fontsize=8)
                    
                    fit_count += 1
                else:
                    print(f"Failed to evaluate model {key}")
                    
            except Exception as e:
                print(f"Error plotting {key}: {e}")

        self.plot_canvas.draw()

    def prev_curve(self):
        if self.current_idx > 0:
            self.load_curve(self.current_idx - 1)

    def next_curve(self):
        if self.current_idx < len(self.curve_list) - 1:
            self.load_curve(self.current_idx + 1)

if __name__ == "__main__":
    root = tk.Tk()
    app = FriedPotatoViewer(root)
    root.mainloop()