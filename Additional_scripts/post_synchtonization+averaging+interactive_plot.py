# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 18:39:53 2021

@author: lpeka
"""

## packages to import

import pandas as pd
import matplotlib.pyplot as plt
import os
import scipy.optimize
import glob
import time
from tkinter import filedialog
import plotly.graph_objects as go
from pastamarkers import markers

global D_list_b
global D_list_o

color_scheme_Solarized_colors = ['#738A05', '#259286', '#2176C7', '#595AB7', '#C61C6F', '#D11C24', '#BD3613', '#A57706']
color_scheme_Solarized_base = ["EAE3CB","FCF4DC","819090","708284","536870","475B62","0A2933","042029"]

""" comparison range  """

F_min=30
F_max=40
interval=1

F_list=list(range(F_min, F_max, interval))
# a list of values according to which the curves are shifted

D_min=900
D_max=2400


# 900 1300 for EMCV
# 1000-1400
interval=1

D_list=list(range(D_min, D_max, interval))


"""functions to be used """



def get_path():
    import_file_path = filedialog.askdirectory()
    return import_file_path
## load a csv file
def getCSV (import_file_path):
    global df
    global F
    global PD_nm
    global filename
    
    df = pd.read_csv (import_file_path)
    #print (df)
    D=df.iloc[:,1] # to select only certain part of the imported table [Row,column]
    F=df.iloc[:,0]
    D_nm=D
    filename=os.path.basename(import_file_path)
    
    return F, D_nm, filename
 
def find_closest(what, where):  
    closest_values=[]
    closest_values_position=[]
    try:
        for i in what:
            absolute_difference_function = lambda cPD : abs(cPD- i)
            value=min(where, key=absolute_difference_function)
            closest_values.append(value)
            where_list=list(where)
            closest_values_position.append(where_list.index(value))
    except:
        pass
    return closest_values_position

def dif_sum(x):
    dif=0
    for i in range(0, len(D_list_b)):
        dif=dif+abs(D_list_b[i]-D_list_o[i]+x)
    
    return dif

def export_data(Force, Distance, name):


    model_data=pd.DataFrame(list(zip(Force, Distance)), columns=['Force [pN]','Distance [nm]'])
    name_model=name[:-4]+'_shifted.csv'
    model_data.to_csv(name_model, index=False, header=True)
    



def export_data_average(Force, Distance, name):


    model_data=pd.DataFrame(list(zip(Force, Distance)), columns=['Force [pN]','Distance [nm]'])
    name_model=name
    model_data.to_csv(name_model, index=False, header=True)    
    
    
    
    
    
    
""" load the curves """




## path to the curves to be shifted
other_path=get_path()
folder_path = other_path + "/*.csv"
print(other_path)
timestamp = time.strftime("%Y%m%d-%H%M%S")
analysis_folder = other_path + '/Postsynchronized_'+ str(F_min) + '-'+str(F_max) +'_' + timestamp
os.mkdir(analysis_folder)

Files = glob.glob(folder_path)


## select the last curve to be the one to which it is alligned

F_b, D_b, name_b = getCSV(Files[-1])
F_list_positions_b=find_closest(F_list, F_b)

D_list_b=[]
for i in F_list_positions_b:
    D_list_b.append(D_b[i])

i=0
while i < len(Files):
    
   
    F_o, D_o, name_o = getCSV(Files[i])
    F_list_positions_o=find_closest(F_list, F_o)
    D_list_o=[]

    for n in F_list_positions_o:
        D_list_o.append(D_o[n])


    

    """ optimize the shift of the other curves """


    ## find optimal shift 

    optim_results=scipy.optimize.minimize(dif_sum, 0)

    x_optim=optim_results['x'][0]



## shift the other curve 

    D_o_new=[]

    for j in D_o:
        D_o_new.append(j-x_optim)


    ## export shifted data 

    export_path=analysis_folder+"\\"+name_o
    print(export_path)
    export_data(F_o, D_o_new, export_path)
    i=i+1

"""start of the averaging part"""




folder_path = analysis_folder  + "/*.csv"

timestamp = time.strftime("%Y%m%d-%H%M%S")
averaged_folder = analysis_folder + '/Averaged_' + timestamp
os.mkdir(averaged_folder)
i=0
Files = glob.glob(folder_path)

D_list_a=[]
F_list_a=[]

for n in range(0, len(D_list)):
    D_list_a.append(0)
    F_list_a.append(0)


while i < len(Files):
    
   
    F_a, D_a, name_a = getCSV(Files[i])
    F_list_positions_a=find_closest(D_list, D_a)

    for j in range(0,len(D_list_a)):
        D_list_a[j]=D_list_a[j]+D_a[F_list_positions_a[j]]
        F_list_a[j]=F_list_a[j]+F_a[F_list_positions_a[j]]
    plt.plot(D_a, F_a)
    i=i+1


F_average=[]
D_average=[]

for k in range(0,len(D_list_a)):
    F_average.append(F_list_a[k]/len(Files))
    D_average.append(D_list_a[k]/len(Files))
    

#plt.plot(D_b, F_b, 'k')
#plt.plot(D_o, F_o, 'r')
    
export_path=averaged_folder+"/averaged.csv"
print(export_path)
export_data_average(F_average, D_average, export_path)


# Initialize an empty list to store dataframes
df_list = []

# Get a list of files in the folder
file_list = os.listdir(analysis_folder)

# Loop through each file in the folder
for file_name in file_list:
    if file_name.endswith(".csv"):  # Check if file is a CSV file
        # Load the CSV file into a dataframe
        file_path = os.path.join(analysis_folder, file_name)
        df = pd.read_csv(file_path)
        
        # Add file name as a column
        df['File'] = os.path.splitext(file_name)[0]
        
        # Append to df_list
        df_list.append(df)

# Concatenate all dataframes in df_list into a single dataframe
combined_df = pd.concat(df_list, ignore_index=True)

# Plot using plotly
fig = go.Figure()

# Counter for indexing color_scheme_Solarized
color_index = 0

for file_name, group in combined_df.groupby("File"):
    fig.add_trace(go.Scatter(x=group["Distance [nm]"], 
                             y=group["Force [pN]"], 
                             mode='lines', 
                             name=file_name, 
                             hoverinfo="name", 
                             showlegend=True,
                             line=dict(color=color_scheme_Solarized_colors[color_index])))
    color_index = (color_index + 1) % len(color_scheme_Solarized_colors)  # Move to the next color in the list

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
fig.write_html(os.path.join(analysis_folder,str(os.path.basename(os.path.normpath(analysis_folder)))+"_force_distance_plot.html"))

# Output file path
print(analysis_folder)
print("Plot saved as force_distance_plot.html")



