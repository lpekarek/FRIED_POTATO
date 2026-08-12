# -*- coding: utf-8 -*-
"""
Lumo-enhanced H5 file viewer with full HF data export
Handles complex BlueLake/H5 structures with deep metadata preservation.
"""

from tkinter import filedialog, Tk, Button, Frame, Label, messagebox, Toplevel
import lumicks.pylake as lk
import os
import matplotlib.pyplot as plt
import h5py
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import SpanSelector
import tempfile
import shutil


class H5DataViewer:
    def __init__(self):
        self.root = Tk()
        self.root.title("H5 Data Viewer - Full Export Mode")
        self.root.geometry("1200x800")
        
        # Data storage
        self.file_new = None
        self.file_h5py = None
        self.filename_new = None
        self.plot_dict = {}       # Stores data for plotting (Force 2x, Distance)
        self.all_hf_data = {}     # Stores ALL HF datasets for export
        self.Force_HF = None
        self.Distance_HF = None
        self.HF_time_ms = None
        self.span_selector = None
        self.selected_region = None
        
        self.setup_gui()
        
    def setup_gui(self):
        """Create the main GUI layout"""
        top_frame = Frame(self.root)
        top_frame.pack(fill='x', padx=10, pady=5)
        
        self.load_btn = Button(top_frame, text="Load H5 File", command=self.load_file, 
                               width=15, height=2)
        self.load_btn.pack(side='left', padx=5)
        
        self.export_btn = Button(top_frame, text="Export Selected Region", 
                                 command=self.export_selected_region, width=15, height=2,
                                 state='disabled')
        self.export_btn.pack(side='left', padx=5)
        
        self.clear_btn = Button(top_frame, text="Clear Selection", 
                                command=self.clear_selection, width=15, height=2)
        self.clear_btn.pack(side='left', padx=5)
        
        self.status_label = Label(top_frame, text="Ready - Load a file to begin", fg="gray")
        self.status_label.pack(side='right', padx=10)
        
        plot_frame = Frame(self.root)
        plot_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        self.fig, (self.ax_force, self.ax_distance) = plt.subplots(2, 1, figsize=(12, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        
        toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        toolbar.update()
        
        instr_frame = Frame(self.root)
        instr_frame.pack(fill='x', padx=10, pady=5)
        instr_label = Label(instr_frame, 
                           text="Instructions: Load File -> Drag on Force plot to select -> Export (slices ALL HF data)",
                           fg="blue")
        instr_label.pack()
        
    def load_file(self):
        """Open file dialog and load ALL HF data"""
        filename = filedialog.askopenfilename(
            title="Select H5 File",
            filetypes=[("H5 files", "*.h5"), ("All files", "*.*")]
        )
        
        if not filename:
            return
            
        try:
            self.status_label.config(text=f"Loading: {os.path.basename(filename)}...")
            self.root.update()
            
            if self.file_h5py:
                self.file_h5py.close()
                
            self.filename_new = filename
            self.file_new = lk.File(filename)
            self.file_h5py = h5py.File(filename, "r")
            
            # Reset storage
            self.all_hf_data = {}
            self.plot_dict = {}
            
            # 1. Load specific data for PLOTTING (Force 2x and Distance)
            if "Force HF/Force 2x" in self.file_h5py:
                self.plot_dict["Force HF/Force 2x"] = self.file_h5py["Force HF/Force 2x"][:]
            if "Distance/Piezo Distance" in self.file_h5py:
                self.plot_dict["Distance/Piezo Distance"] = self.file_h5py["Distance/Piezo Distance"][:]
            
            if not self.plot_dict:
                raise ValueError("Could not find Force HF/Force 2x or Distance/Piezo Distance")
                
            self.Force_HF = self.plot_dict["Force HF/Force 2x"]
            self.Distance_HF = self.plot_dict["Distance/Piezo Distance"]
            
            # Calculate time array based on the main Force HF data
            HF_time_ns = np.arange(0, len(self.Force_HF))
            self.HF_time_ms = HF_time_ns / 10**6
            
            # 2. Load ALL HF datasets for EXPORT
            target_length = len(self.Force_HF)
            
            def collect_hf_data(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if (len(obj.shape) == 1 and 
                        obj.dtype in [np.float64, np.float32] and 
                        obj.size == target_length):
                        if ("Force HF" in name or "Trap position" in name or 
                            "Piezo" in name or "Corrected" in name):
                            self.all_hf_data[name] = obj[:]
            
            self.file_h5py.visititems(collect_hf_data)
            
            self.status_label.config(
                text=f"Loaded: {len(self.Force_HF)} samples | Found {len(self.all_hf_data)} HF datasets"
            )
            print(f"Loaded HF datasets for export: {list(self.all_hf_data.keys())}")
            
            # Call the METHOD (renamed to avoid collision)
            self.plot_data()
            self.export_btn.config(state='normal')
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
            self.status_label.config(text="Error loading file")
            import traceback
            print(traceback.format_exc())
            
    def plot_data(self):
        """Plot the main HF data"""
        self.ax_force.clear()
        self.ax_distance.clear()
        
        # Plot Force 2x
        if "Force HF/Force 2x" in self.plot_dict:
            self.ax_force.plot(self.HF_time_ms, self.plot_dict["Force HF/Force 2x"], c='b', 
                              linewidth=0.5, markersize=0.5, marker="o", label='Force 2x')
            self.ax_force.set_ylabel("Force, pN")
            self.ax_force.legend(loc='upper left')
            self.ax_force.grid(True, alpha=0.3)
        
        # Plot Distance
        if "Distance/Piezo Distance" in self.plot_dict:
            self.ax_distance.plot(self.HF_time_ms, self.plot_dict["Distance/Piezo Distance"] * 1000, 
                                 c='b', linewidth=0.5, markersize=0.5, marker="o", label='Distance')
            self.ax_distance.set_ylabel("Distance, nm")
            self.ax_distance.set_xlabel("Time, ms")
            self.ax_distance.legend(loc='upper left')
            self.ax_distance.grid(True, alpha=0.3)
        
        self.fig.suptitle(f"H5 Data: {os.path.basename(self.filename_new)}", fontsize=12)
        self.fig.tight_layout()
        
        self.add_span_selector()
        self.canvas.draw()
        
    def add_span_selector(self):
        """Add interactive span selector"""
        if self.span_selector:
            self.span_selector.disconnect()
            
        def onselect(xmin, xmax):
            idx_min = np.searchsorted(self.HF_time_ms, xmin)
            idx_max = np.searchsorted(self.HF_time_ms, xmax)
            
            self.selected_region = {
                'xmin': xmin,
                'xmax': xmax,
                'idx_min': idx_min,
                'idx_max': idx_max,
                'duration_ms': xmax - xmin
            }
            
            self.status_label.config(
                text=f"Selected: {xmin:.2f}ms - {xmax:.2f}ms ({self.selected_region['duration_ms']:.2f}ms)"
            )
            
        self.span_selector = SpanSelector(
            self.ax_force,
            onselect,
            "horizontal",
            useblit=True,
            props=dict(alpha=0.5, facecolor="red"),
            interactive=True,
            drag_from_anywhere=True
        )
        
    def clear_selection(self):
        self.selected_region = None
        self.status_label.config(text="Selection cleared")
        
    def export_selected_region(self):
        """Export selected region, slicing ALL HF datasets and preserving ALL metadata"""
        if not self.selected_region:
            messagebox.showwarning("Warning", "Please select a region first")
            return
            
        idx_min = self.selected_region['idx_min']
        idx_max = self.selected_region['idx_max']
        duration_samples = idx_max - idx_min
        
        output_filename = filedialog.asksaveasfilename(
            title="Save Exported Region",
            defaultextension=".h5",
            filetypes=[("H5 files", "*.h5")]
        )
        
        if not output_filename:
            return
            
        try:
            fd, temp_path = tempfile.mkstemp(suffix='.h5')
            os.close(fd)
            
            try:
                # STEP 1: Deep copy the ENTIRE original file structure
                with h5py.File(self.filename_new, 'r') as src, h5py.File(temp_path, 'w') as dst:
                    self.deep_copy_with_attrs(src, dst)
                
                # STEP 2: Open temp file and slice the HF datasets
                with h5py.File(temp_path, 'r+') as f:
                    sliced_count = 0
                    
                    for ds_path, original_data in self.all_hf_data.items():
                        if ds_path not in f:
                            print(f"Warning: Path {ds_path} not found in copied file structure")
                            continue
                        
                        try:
                            orig_ds = f[ds_path]
                            
                            if orig_ds.size != len(original_data):
                                print(f"Skipping {ds_path}: Size mismatch ({orig_ds.size} vs {len(original_data)})")
                                continue
                                
                            # Get properties before deleting
                            dtype = orig_ds.dtype
                            chunks = orig_ds.chunks if orig_ds.chunks else None
                            compression = orig_ds.compression
                            compression_opts = orig_ds.compression_opts
                            
                            # Slice the data
                            new_data = original_data[idx_min:idx_max]
                            
                            # Delete old and create new
                            del f[ds_path]
                            f.create_dataset(ds_path, data=new_data, 
                                           dtype=dtype, chunks=chunks,
                                           compression=compression, 
                                           compression_opts=compression_opts)
                            
                            # Re-apply attributes from the original source file
                            with h5py.File(self.filename_new, 'r') as src_check:
                                if ds_path in src_check:
                                    src_ds = src_check[ds_path]
                                    for attr_name in src_ds.attrs.keys():
                                        try:
                                            f[ds_path].attrs[attr_name] = src_ds.attrs[attr_name]
                                        except Exception as e:
                                            print(f"Warning: Could not copy attr '{attr_name}' to {ds_path}: {e}")
                            
                            sliced_count += 1
                            print(f"Sliced: {ds_path} ({len(new_data)} samples)")
                            
                        except Exception as e:
                            print(f"Error slicing {ds_path}: {e}")
                    
                    # Update global metadata if needed
                    if 'num_samples' in f.attrs:
                        f.attrs['num_samples'] = duration_samples
                    
                    print(f"Successfully sliced {sliced_count} datasets.")
                    
                # Move to final destination
                shutil.move(temp_path, output_filename)
                
                self.status_label.config(text=f"Exported to: {os.path.basename(output_filename)}")
                messagebox.showinfo("Success", 
                    f"Exported {duration_samples} samples from {sliced_count} datasets!\nMetadata preserved.")
                
            finally:
                if os.path.exists(temp_path):
                    try:
                        os.remove(temp_path)
                    except:
                        pass

        except Exception as e:
            messagebox.showerror("Error", f"Failed to export:\n{str(e)}")
            import traceback
            print(traceback.format_exc())

    def deep_copy_with_attrs(self, source, dest):
        """Recursively copy groups, datasets, and ALL attributes with safe chunking"""
        def _copy_recursive(source_item, dest_parent, name):
            if isinstance(source_item, h5py.Group):
                if name not in dest_parent:
                    dest_group = dest_parent.create_group(name)
                else:
                    dest_group = dest_parent[name]
                
                # Copy group attributes
                for attr_name, attr_val in source_item.attrs.items():
                    try:
                        dest_group.attrs[attr_name] = attr_val
                    except Exception as e:
                        print(f"Warning: Failed to copy group attr '{attr_name}': {e}")
                
                for key in source_item.keys():
                    _copy_recursive(source_item[key], dest_group, key)
                    
            elif isinstance(source_item, h5py.Dataset):
                # Handle Scalars
                if source_item.shape == ():
                    if name not in dest_parent:
                        dest_parent.create_dataset(name, data=source_item[()], dtype=source_item.dtype)
                    # Copy attributes
                    for attr_name, attr_val in source_item.attrs.items():
                        try:
                            dest_parent[name].attrs[attr_name] = attr_val
                        except Exception as e:
                            print(f"Warning: Failed to copy dataset attr '{attr_name}': {e}")
                    return

                # Handle Arrays
                if name not in dest_parent:
                    dtype = source_item.dtype
                    shape = source_item.shape
                    original_chunks = source_item.chunks
                    
                    # --- FIX: Adjust Chunk Size ---
                    new_chunks = None
                    if original_chunks:
                        # If original chunks are larger than the data shape, reduce them
                        # We take the minimum of (original_chunk_dim, data_dim)
                        # But ensure we don't create chunks of size 0 or 1 if possible (keep reasonable defaults)
                        adjusted_chunks = tuple(
                            min(c, s) if s > 0 else c 
                            for c, s in zip(original_chunks, shape)
                        )
                        
                        # Ensure chunks are not larger than data
                        if all(c <= s for c, s in zip(adjusted_chunks, shape)):
                            new_chunks = adjusted_chunks
                        else:
                            # Fallback: create a single chunk of the full data size
                            new_chunks = shape
                    
                    # Get compression settings
                    compression = source_item.compression
                    compression_opts = source_item.compression_opts
                    shuffle = source_item.shuffle if hasattr(source_item, 'shuffle') else False
                    fletcher32 = source_item.fletcher32 if hasattr(source_item, 'fletcher32') else False
                    maxshape = source_item.maxshape if hasattr(source_item, 'maxshape') else None

                    try:
                        dest_parent.create_dataset(
                            name, 
                            data=source_item[:], 
                            dtype=dtype, 
                            shape=shape,
                            chunks=new_chunks,
                            compression=compression,
                            compression_opts=compression_opts,
                            shuffle=shuffle,
                            fletcher32=fletcher32,
                            maxshape=maxshape
                        )
                    except Exception as e:
                        # If specific chunking fails, try without chunks (contiguous)
                        print(f"Warning: Chunking failed for {name}, trying contiguous: {e}")
                        dest_parent.create_dataset(
                            name, 
                            data=source_item[:], 
                            dtype=dtype, 
                            shape=shape,
                            compression=compression,
                            compression_opts=compression_opts
                        )
                
                # Copy dataset attributes
                for attr_name, attr_val in source_item.attrs.items():
                    try:
                        dest_parent[name].attrs[attr_name] = attr_val
                    except Exception as e:
                        print(f"Warning: Failed to copy dataset attr '{attr_name}': {e}")

        # Copy root attributes
        for attr_name, attr_val in source.attrs.items():
            try:
                dest.attrs[attr_name] = attr_val
            except Exception as e:
                print(f"Warning: Failed to copy root attr '{attr_name}': {e}")
        
        for key in source.keys():
            _copy_recursive(source[key], dest, key)

    def run(self):
        self.root.mainloop()
        
    def cleanup(self):
        if self.file_h5py:
            self.file_h5py.close()
        self.root.destroy()


if __name__ == "__main__":
    viewer = H5DataViewer()
    try:
        viewer.run()
    finally:
        viewer.cleanup()