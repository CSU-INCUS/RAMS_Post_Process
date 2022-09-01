# RAMS_Post_Process
Summary: Code to calculate variables and figures from the RAMS model files that are necessary for the INCUS forward model

This code is comprised of a driver file (postprocess_rams.py) and a file with functions used by the driver (fx_postprocess_rams.py).
The user needs to specify 4 paths/filenames in the driver file (the RAMS datafile, the RAMS headerfile, and the path for where to save the new data and path to save figures).
The code has been developed for the NASA INvestigation of Convective UpdraftS (INCUS) mission, and its goal is to provide necessary information for the INCUS forward model.
Specifically, this code calculates 2-dimensional (x,y) variables of 10m wind speed and surface temperature from the native RAMS variables and saves these varibles to a new netcdf file using xarray.
This code also create a summary figure of the integrated total condensate and vertical velocities from the RAMS data file for the entire domain.
