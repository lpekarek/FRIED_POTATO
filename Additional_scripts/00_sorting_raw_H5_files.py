import tkinter as tk
from tkinter import filedialog
import glob
import os
from pathlib import Path
import time

def getRAW_folder():
    global Files
    global filename_i
    global folder_name
    global forward_folder
    global reverse_folder
    global folder

    folder = filedialog.askdirectory()
    folder_path = str(folder + "/*.h5")
    Files = glob.glob(folder_path)
    print('Files to analyse:', len(Files))

    if len(Files) == 0:
        print('No file of the selected data type in the folder!')
        return

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    print("Timestamp =", timestamp)
    print('Start of analysis:', timestamp)
    
    # Save the analysed data to a newly created folder with a timestamp
    folder_name = os.path.basename(folder)
    forward_folder = os.path.join(folder, f'{folder_name}_FD_forward')
    reverse_folder = os.path.join(folder, f'{folder_name}_FD_reverse')
    scan_folder = os.path.join(folder, f'{folder_name}_scans')
    kymograph_folder = os.path.join(folder, f'{folder_name}_kymographs')
    marker_folder = os.path.join(folder, f'{folder_name}_markers')

    for i in [forward_folder,reverse_folder,scan_folder , kymograph_folder, marker_folder]:
        try:
            os.mkdir(i)
        
        except:
            pass

    for file_path in Files:
        filename_i = os.path.basename(file_path)
        print(f'Processing file: {filename_i}')

        if "Scan" in filename_i:
            target_path = os.path.join(scan_folder, filename_i)
            print(f'Moving {filename_i} to {target_path}')
            Path(file_path).rename(target_path)
        elif "Kymograph" in filename_i:
            target_path = os.path.join(kymograph_folder, filename_i)
            print(f'Moving {filename_i} to {target_path}')
            Path(file_path).rename(target_path)
        elif "Marker" in filename_i:
            target_path = os.path.join(marker_folder, filename_i)
            print(f'Moving {filename_i} to {target_path}')
            Path(file_path).rename(target_path)
        elif "FD Curve" in filename_i and "r.h5" in filename_i:
            target_path = os.path.join(reverse_folder, filename_i)
            print(f'Moving {filename_i} to {target_path}')
            Path(file_path).rename(target_path)
        elif "FD Curve" in filename_i and "f.h5" in filename_i:
            target_path = os.path.join(forward_folder, filename_i)
            print(f'Moving {filename_i} to {target_path}')
            Path(file_path).rename(target_path)
        else:
            print(f'File {filename_i} does not match any criteria.')

if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    getRAW_folder()
