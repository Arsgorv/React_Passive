#!/usr/bin/env python3
"""
Converts Open Ephys event .npy files in a folder to a single .mat file.
Run as:
    python convertEvents2Mat.py -p /path/to/session/events
"""

import os
import sys
import argparse

# --- Dependency check ---
try:
    import numpy as np
    from scipy.io import savemat
except ImportError as e:
    print("Missing required packages. Please run:")
    print("    pip install numpy scipy")
    sys.exit(1)

# --- Argument parser ---
parser = argparse.ArgumentParser(description="Converts open-ephys events to mat files")
parser.add_argument('-p','--path', type=str, required=True,
                    help="Path to 'events' folder of recording session")
args = parser.parse_args()

path = os.path.abspath(os.path.expanduser(args.path))

# --- Path validity check (optional, can be commented out if not needed) ---
if '.' in path:
    if 'Rhythm_FPGA-100.0' not in path:
        raise ValueError('Please give full path, not relative path.')

def findFolderName(path):
    """Returns last folder in path"""
    return os.path.basename(os.path.normpath(path))

def findPreviousFolderName(path):
    """Returns parent folder name"""
    parent = os.path.dirname(os.path.normpath(path))
    return os.path.basename(parent)

# --- Main conversion logic ---
saved_files = 0

for r, d, f in os.walk(path):
    # Only process folders containing .npy files
    npy_files = [file for file in f if file.endswith('.npy')]
    if len(npy_files) == 0:
        continue

    data = {}
    for file in npy_files:
        file_path = os.path.join(r, file)
        try:
            arr = np.load(file_path)
            key = os.path.splitext(file)[0]
            data[key] = arr
        except Exception as e:
            print(f"Error loading {file_path}: {e}")

    if data:
        # Save .mat file to same folder as the .npy files, named ParentFolder_CurrentFolder.mat
        current_folder = findFolderName(r)
        parent_folder = findPreviousFolderName(r)
        mat_filename = f"{parent_folder}_{current_folder}.mat"
        mat_filepath = os.path.join(path, mat_filename)

        try:
            savemat(mat_filepath, data)
            print(f"Saved {mat_filepath} (contains: {list(data.keys())})")
            saved_files += 1
        except Exception as e:
            print(f"Error saving {mat_filepath}: {e}")

if saved_files == 0:
    print("No .npy files found for conversion in the specified path.")
else:
    print(f"Conversion complete. {saved_files} .mat file(s) written.")

