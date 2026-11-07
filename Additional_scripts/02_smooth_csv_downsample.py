# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 16:13:10 2021

@author: lpeka
"""
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import filedialog
from tkinter import ttk
import pandas as pd
import time
import glob
import os
import numpy as np
from pathlib import Path


DS_factor=100





def getRAW_folder():
    global Files
    global filename_i
    global directory_i
    global analysis_folder
    global timestamp
    global figure1
    global x_STD_1
    global F_min
    global Force_Distance
    global filteredDistance_ready
    global filteredForce

    folder = filedialog.askdirectory()
    folder_path = str(folder + "/*.csv")
    #print(folder_path)
    Files = glob.glob(folder_path)
    print('files to analyse', len(Files))

    print('Files to analyse: ' + str(len(Files)))
    

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    print("Timestamp = " + timestamp)
    print('Start of analysis: ' + str(timestamp))
    

    if len(Files) == 0:
        print('No file of the selected data type in the folder!')
       
    else:
        print('Please do not close the program!')

    #save the analysed data to a new created folder with a timestamp
    analysis_folder = str(folder+'/Downsampled_by_'+ str(DS_factor) + '_'+ timestamp)
    os.mkdir(analysis_folder)

    #iterate through the h5 files in the folder
    i = 0
    total_results=[]


    while i < len(Files):
        print(Files[i])
        if "smooth" in Files[i]: #print(Files[i]):
            directory_i = Path(Files[i])
            filename_i = directory_i.name[:-4]
            downsampleCSV(Files[i])
            #print(df)
        i=i+1


        
        
        


def downsampleCSV (import_file_path):
    global df
    global F
    global PD_nm
    global filename_i
    global DS_factor
    global F_ds
    
    df = pd.read_csv (import_file_path)
    #print (df)
    PD=df.iloc[:,1] # to select only certain part of the imported table [Row,column]
    F=df.iloc[:,0]
    PD_nm=PD*1000
    
    F_ds= []
    PD_ds=[]
    n=0
    while n < len(F):
        F_ds.append(F[n])
        PD_ds.append(PD_nm[n])
        n=n+DS_factor
        
    print(F_ds)
    export_data(PD_ds, F_ds)
        
def export_data(PD_ds, F_ds):    
    global name_DS_data
    DS_data=pd.DataFrame(list(zip(F_ds, PD_ds)), columns=[ 'Force [pN]','Distance [nm]'])

    name_DS_data=analysis_folder + "/" + filename_i+'_downsamplesd.csv'
    DS_data.to_csv(name_DS_data, index=False, header=True)
        

root = tk.Tk()
    
getRAW_folder()  


  
    