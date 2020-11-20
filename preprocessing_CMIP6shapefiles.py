#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 12:01:04 2020

@author: danielfurman
"""
# This file unzips climate tifs downloaded from Worldclim.org, and
# the expanded files are saved to a HardDrive folder "output_dir".
# Manually unzipping eight tif files for a single ssp scenario
# took 3 minutes and 24 seconds, in part due to an extensive opening
# of subfolders. In contrast, this code runs in approximately 31
# seconds, speeding processing by over 600%. Overall, we saved over
# 30 minutes when we downloaded CMIP6 data for 12 scenarios.

# see https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html

# It is important the files match the eight CMIP6 models in
# model_names, and that the data are 2.5m resolution.

import time  # we will calculate the runtime
import zipfile  # library for unzipping function
import glob  # data loading functions
import shutil  # moving files within directories functions

# output folder
ssp = 'ssp126'
output_dir = '/volumes/HardDrive/' + ssp + '_2061-2080'

# make sure only 8 relevant zip are in downloads
filenames = sorted(glob.glob('/Users/danielfurman/Downloads/wc2.1*.zip'))
print('The file names are:', filenames)  # print filenames
len(filenames)  # Confirmed 8 total files

start_time = time.time()

for f in filenames:  # f is index with string filenames
    with zipfile.ZipFile(f, "r") as zip_ref:
        zip_ref.extractall(output_dir)

model_names = ['BCC-CSM2-MR', 'CanESM5', 'CNRM-CM6-1', 'CNRM-ESM2-1',
               'IPSL-CM6A-LR', 'MIROC-ES2L', 'MIROC6', 'MRI-ESM2-0']

for model in model_names:
    original = glob.glob(output_dir +
                         '/share/spatial03/worldclim/cmip6/7_fut/2.5m/' +
                         model + '/' + ssp + '/wc*.tif')
    shutil.move(original[0], output_dir)

print("--- %s seconds ---" % (time.time() - start_time))
