

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 17:27:17 2024

@author: lupe184g
"""
import os
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import plotly.graph_objects as go

"""_Solarized_"""
color_scheme_colors = ['#738A05', '#259286', '#2176C7', '#595AB7', '#C61C6F', '#D11C24', '#BD3613', '#A57706']
color_scheme_base = ["#EAE3CB","#FCF4DC","#819090","#708284","#536870","#475B62","#0A2933","#042029"]
"""_Nord_"""
#color_scheme_colors=["#2A6366","#B48EAD","#A3BE8C","#EBCB8B","#D08770","#BF616A","#5E81AC","#81A1C1","#88C0D0","#8FBCBB"]
#color_scheme_base=["#ECEFF4","#E5E9F0","#D8DEE9","#4C566A","#434C5E","#3B4252","#2E3440"]


root = tk.Tk()

def openFile():
    global folder_path_f
    global folder_path_r
    # Specify the folder path containing CSV files
    folder_path_f = filedialog.askdirectory()
    folder_path_r = filedialog.askdirectory()
    root.destroy()

openFile()

# Initialize an empty list to store dataframes
df_list_f = []
df_list_r = []
# Get a list of files in the folder
file_list_f = os.listdir(folder_path_f)
file_list_r = os.listdir(folder_path_r)

# Loop through each file in the folder
for file_name in file_list_f:
    if file_name.endswith(".csv"):  # Check if file is a CSV file
        # Load the CSV file into a dataframe
        file_path = os.path.join(folder_path_f, file_name)
        df = pd.read_csv(file_path)
        
        # Add file name as a column
        df['File'] = os.path.splitext(file_name)[0]
        
        # Append to df_list
        df_list_f.append(df)

# Loop through each file in the folder
for file_name in file_list_r:
    if file_name.endswith(".csv"):  # Check if file is a CSV file
        # Load the CSV file into a dataframe
        file_path = os.path.join(folder_path_r, file_name)
        df = pd.read_csv(file_path)
        
        # Add file name as a column
        df['File'] = os.path.splitext(file_name)[0]
        
        # Append to df_list
        df_list_r.append(df)

# Concatenate all dataframes in df_list into a single dataframe
combined_df_f = pd.concat(df_list_f, ignore_index=True)
combined_df_r = pd.concat(df_list_r, ignore_index=True)
# Plot using plotly
fig = go.Figure()

# Counter for indexing color_scheme_Solarized
color_index = 0

for file_name, group in combined_df_f.groupby("File"):
    fig.add_trace(go.Scatter(x=group["Distance [nm]"], 
                             y=group["Force [pN]"], 
                             mode='lines', 
                             name=file_name, 
                             hoverinfo="name", 
                             showlegend=True,
                             line=dict(color=color_scheme_colors[color_index])))
    color_index = (color_index + 1) % len(color_scheme_colors)  # Move to the next color in the list

color_index = 0
for file_name, group in combined_df_r.groupby("File"):
    fig.add_trace(go.Scatter(x=group["Distance [nm]"], 
                             y=group["Force [pN]"], 
                             mode='lines', 
                             name=file_name, 
                             hoverinfo="name", 
                             showlegend=True,
                             line=dict(color=color_scheme_base[color_index])))
    color_index = (color_index + 1) % len(color_scheme_base)  # Move to the next color in the list









fig.update_layout(title="Force vs Distance", 
                  xaxis_title="Distance, nm", 
                  hoverlabel_namelength=-1,
                  yaxis_title="Force, pN", 
                  template="plotly_white", 
                  legend=dict(traceorder="normal"))

# Customize hover effects
fig.update_traces(hovertemplate=None)

# Check if the plot is generated
print("Plot information:")
print(fig)

# Save the plot to an HTML file
fig.write_html(os.path.join(folder_path_f, "force_distance_plot_f_and_r.html"))

# Output file path
print("Plot saved as " + str(os.path.join(folder_path_f, "force_distance_plot_f_and_r.html")))