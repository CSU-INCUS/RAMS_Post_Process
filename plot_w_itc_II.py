#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: petermarinescu

Purpose: to create quick looks of the 30s, LES output

Notes: Code is parallized with jug; Need to update path and savepath varibles below
       This script is run with the run_plot_w_itc_II.sh and ple_run_jug_plot_w_itc scripts 

"""

# Import functions from fx_postproc_RAMS.py file in the same directory
from fx_postproc_RAMS import read_head, calc_surf_temp, calc_10m_wind, plot_w_itc, save_rams2D_to_netcdf
import glob
import numpy as np
from jug import TaskGenerator

path = '/nobackup/pmarines/PROD/ARG1.2-R/G3/out_30s/'
savefig_path = '/nobackup/pmarines/PROD/ARG1.2-R/G3/out_30s/Plots/'
run_ID = 'ARG1.2-R-V1-G3'

# Specify file suffix and use glob to get list of files
h5filepath = path+'a-L*g3.h5'
files = sorted(glob.glob(h5filepath))
hefilepath = path+'a-L*head.txt'
hfiles = sorted(glob.glob(hefilepath))

# Check to make sure there is a header file for each RAMS hdf file
if len(files) != len(hfiles):
    print('Not equal header and data files')
    exit()

# Read RAMS header file to get coordinate information
zm, zt, nx, ny, dxy, npa = read_head(hfiles[0],files[0])

# Create function loop_files to help with parallization
#@TaskGenerator
def loop_files(files,hfiles,zt,savefig_path):

    for i in np.arange(0,len(files)):
    # Specify filename (including path) of the RAMS datafile and associated RAMS headerfile    

        print(i)
        datafile = files[i]
        headfile = hfiles[i]
        #savedata_path = '/Users/petermarinescu/Research/INCUS/Code/Test_Data/' # Path to save netcdf files with data

        # Read RAMS header file to get coordinate information
        zm, zt, nx, ny, dxy, npa = read_head(headfile,datafile)

        # Call function to calculate surface temperature (proxy for skin temperature)
        # This ultimately takes an area-weighted average of the vegetation temperature and soil temperature
        #temp_surf = calc_surf_temp(datafile)

        # Call function to calculate 10 m wind speed 
        # This is based on Monin–Obukhov theory and follows the  structure of the RAMS revu code
        #wind_10m = calc_10m_wind(datafile,zt)

        # Save these two variables to a netcdf file using xarray at savedata_path
        #save_rams2D_to_netcdf(datafile, ['temp_surf','wind_10m'], [temp_surf,wind_10m], savedata_path)

        # Plot maximum vertical velocity and vertically integrated condensate over entire domain
        plot_w_itc(datafile,savefig_path,zt)

# Run loop files function
loop_files(files,hfiles,zt,savefig_path,run_ID)
