#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: petermarinescu

Purpose: to post-processes RAMS model output to create 
         inputs needed for forward simulator and
         figures to help visualize where updrafts/clouds are
         for each model snapshot time

Notes: Code is primarly based on functions from RAMS-revu code

"""

# Import functions from fx_postproc_RAMS.py file in the same directory
from fx_postproc_RAMS import read_head, calc_surf_temp, calc_10m_wind, plot_w_itc, save_rams2D_to_netcdf

# Specify filename (including path) of the RAMS datafile and associated RAMS headerfile 
datafile = '/Users/petermarinescu/Research/INCUS/Code/Test_ACPC/a-A-2013-06-19-160000-g3.h5'
headfile = '/Users/petermarinescu/Research/INCUS/Code/Test_ACPC/a-A-2013-06-19-160000-head.txt'
savedata_path = '/Users/petermarinescu/Research/INCUS/Code/Test_Data/' # Path to save netcdf files with data
savefig_path = '/Users/petermarinescu/Research/INCUS/Code/Test_Figure/' # Path to save summary figure

# Read RAMS header file to get coordinate information
zm, zt, nx, ny, dxy, npa = read_head(headfile,datafile)

# Call function to calculate surface temperature (proxy for skin temperature)
# This ultimately takes an area-weighted average of the vegetation temperature and soil temperature
temp_surf = calc_surf_temp(datafile)

# Call function to calculate 10 m wind speed 
# This is based on Monin–Obukhov theory and follows the  structure of the RAMS revu code
wind_10m = calc_10m_wind(datafile,zt)

# Save these two variables to a netcdf file using xarray at savedata_path
save_rams2D_to_netcdf(datafile, ['temp_surf','wind_10m'], [temp_surf,wind_10m], savedata_path)

# Plot maximum vertical velocity and vertically integrated condensate over entire domain
plot_w_itc(datafile,savefig_path,zt)