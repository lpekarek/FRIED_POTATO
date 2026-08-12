# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:27:37 2024

@author: lupe184g
"""


from tkinter import filedialog
import lumicks.pylake as lk
import pandas as pd
import os
import matplotlib.pyplot as plt
import h5py
import numpy as np

"""loading files """

#filename_new="Z:\\OT_RNA_folding\\OT_RAW\\2025-06-11\\T01_M1\\20250611-112425 FD Curve HD2_B04_P000_R01_T01_Multi.h5"
filename_new = "C:\\Users\\lupe184g\\Scripts_LP\\FRIED_POTATO_troubleshooting\\Kelly_data\\20260510-154441 Marker HOP1b-3_f3_2.5kbH_BP1_M (4).h5"
#filename_new = "D:/Kelly/20260510-154441 Marker HOP1b-3_f3_2.5kbH_BP1_M (4).h5"
print('filename')
file_new=lk.File(filename_new)

file_h5py=h5py.File(filename_new, "r")
print('file loaded')
Force_HF=file_h5py["Force HF/Force 2x"][:]
Distance_HF = file_h5py["Distance/Piezo Distance"][:]
print(len(Force_HF))
print(len(Distance_HF))

#Force_LF=file_h5py["Force LF/Force 2x"][:]
#Distance_LF = file_h5py["Distance/Distance 1"][:]
#print(len(Force_LF))
#print(len(Distance_LF))

#time in ms
#LF_time_ms=(file_h5py["Force LF/Force 2x"]["Timestamp"]-file_h5py["Force LF/Force 2x"]["Timestamp"][0])/10**6
#time_length_ns=(file_h5py["Force HF/Force 2x"]["Timestamp"][-1]-file_h5py["Force HF/Force 2x"]["Timestamp"][0])
ns_between_HF_data_taken=1/78100*10**9
#ns_between_HF_data_taken_2=time_length_ns/len(Force_HF)
HF_time_ns=np.arange(0, len(Force_HF))
HF_time_ms=HF_time_ns/10**6


""" ploting the FD data"""



# x = range(100)
# y = range(100,200)
# fig = plt.figure()
# ax1 = fig.add_subplot(111)

# ax1.plot(HF_time_ms, Force_HF, c='b', linewidth=0.5, markersize=0.5, marker="o", label='HF')
# ax1.plot(LF_time_ms,Force_LF['Value'], c='r',linewidth=0.5,markersize=0.5,  marker="o", label='LF')
# plt.legend(loc='upper left')

# plt.xlabel("time, ms")
# plt.ylabel("Force, pN")

# ax2 = fig.add_subplot(222)
# ax2.plot(HF_time_ms, Distance_HF, c='b', linewidth=0.5, markersize=0.5, marker="o", label='HF')
# ax2.plot(LF_time_ms,Distance_LF['Value'], c='r',linewidth=0.5,markersize=0.5,  marker="o", label='LF')

# plt.show()



# Initialise the subplot function using number of rows and columns 
figure, axis = plt.subplots(2)

axis[0].plot(HF_time_ms, Force_HF, c='b', linewidth=0.5, markersize=0.5, marker="o", label='HF') 
#axis[0].plot(LF_time_ms, Force_LF['Value'], c='r',linewidth=1, label='LF')
#axis[0].set_title("Force") 
axis[0].set(ylabel="Force, pN")
axis[0].legend(loc='upper left')

axis[1].plot(HF_time_ms, Distance_HF*1000, c='b', linewidth=0.5, markersize=0.5, marker="o", label='HF')
#axis[1].plot(LF_time_ms,Distance_LF['Value']*1000, c='r',linewidth=1, label='LF')
#axis[1].set_title("Distance") 
axis[1].set(ylabel="Distance, nm")
plt.xlabel("time, ms")

# plt.ylabel("Force, pN")

plt.show()
