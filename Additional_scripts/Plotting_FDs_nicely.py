import pandas as pd
import matplotlib.pyplot as plt
import os
import tkinter as tk
from tkinter import filedialog
import matplotlib.patches as patches

# Define the directories containing the CSV files
#directory1 = r'C:\Users\lupe184g\Desktop\Postdoc\01_Projects\OT_data\2024_05_04_data_for_MJ_plot\Representative_FDs\'
#directory2 = r'C:\Users\lupe184g\Desktop\Postdoc\01_Projects\OT_data\2024_05_04_data_for_MJ_plot\Representative_FDs\forwards_aligned\CSPC'

root = tk.Tk()

def openFile():
    global directory1
    global directory2
    # Specify the folder path containing CSV files
    directory1 = filedialog.askdirectory()
    directory2 = filedialog.askdirectory()
    root.destroy()

openFile()

# Define the specific hex color codes for each folder
color_hex1 = '#CC0078'  # Pink color for the first folder #RBP #3883CB #RNA #CC0078
color_hex2 = '#667C85'  # Blue color for the second folder #grey reverse #667C85 ! #YBX1 #778f3e
font_size = 20
Plot_title=r"HOTAIR D2"  #+- 1 $\mu$M CspC 
data1_name = "unfolding"
data2_name = "refolding"

line_thickness = 2
min_x_value, max_x_value = 800, 1449
min_y_value, max_y_value = -1 , 44

# Add a scale bar
scalebar_length = 100  # Length of the scale bar in x-axis units
scalebar_height = 0.01 * (max_y_value - min_y_value)  # Relative height of the scale bar
scalebar_x_position = max_x_value-scalebar_length*1.001  # Position of the scale bar on the x-axis
scalebar_y_position = min_y_value + scalebar_height * 2  # Position on the y-axis


data_name = Plot_title


# Function to load CSV data
def load_csv(filename):
    return pd.read_csv(filename, usecols=['Force [pN]', 'Distance [nm]'])


# Function to plot data from a specific directory
def plot_data_from_directory(directory, color_hex, label, initial_alpha=1.0, alpha_step=0.15):
    files = sorted(os.listdir(directory))
    for i, file in enumerate(files):
        if file.endswith('.csv'):
            filepath = os.path.join(directory, file)
            data = load_csv(filepath)
            # Add a label only for the first file in the directory
            if i == 0:
                plt.plot(data['Distance [nm]'], data['Force [pN]'], label=label,
                         color=color_hex, alpha=initial_alpha, linewidth=line_thickness)
            else:
                plt.plot(data['Distance [nm]'], data['Force [pN]'],
                         color=color_hex, alpha=initial_alpha - i * alpha_step, linewidth=line_thickness)

# Plotting
plt.figure(figsize=(12, 8), dpi=600)
plot_data_from_directory(directory1, color_hex1, data1_name)
plot_data_from_directory(directory2, color_hex2, data2_name)

plt.title(Plot_title, fontsize=font_size)

# Set the tick and axis properties
ax = plt.gca()  # Get current axis
ax.spines['top'].set_linewidth(2)    # Set the width of the top spine
ax.spines['right'].set_linewidth(2)  # Set the width of the right spine
ax.spines['left'].set_linewidth(2)   # Set the width of the left spine
ax.spines['bottom'].set_linewidth(2) # Set the width of the bottom spine
ax.spines['top'].set_color('white') 
ax.spines['right'].set_color('white')


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

plt.xlabel('Relative distance, nm', fontsize=font_size)
plt.ylabel('Force, pN', fontsize=font_size)
plt.legend(fontsize=font_size)


save_dir = r'C:\Users\lupe184g\Desktop\Postdoc\01_Projects\OT_data\Fiona\2025_01_24_FDs_for_nice_plots'
os.makedirs(save_dir, exist_ok=True)  # This creates the directory if it doesn't exist

file_path_png = os.path.join(save_dir, f'{data_name}_FD_curves_zoomed.png')
file_path_svg = os.path.join(save_dir, f'{data_name}_FD_curves_zoomed.svg')

plt.savefig(file_path_png)  # Save as PNG
plt.savefig(file_path_svg, format='svg')  # Save as SVG


#plt.savefig(r'C:\Users\lupe184g\Desktop\Postdoc\01_Projects\OT_data\2024_05_10_data_for_MJ_plot\Plots\ '+str(data_name)+ '_FD_curves'+'.png')  # Save as PNG
#plt.savefig(r'C:\Users\lupe184g\Desktop\Postdoc\01_Projects\OT_data\2024_05_10_data_for_MJ_plot\Plots\ '+str(data_name)+ '_FD_curves'+'.svg', format='svg')  # Save as svg



#plt.show()
