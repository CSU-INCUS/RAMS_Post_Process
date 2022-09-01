# RAMS_Post_Process
Summary: Code to calculate variables and figures from the RAMS model files that are necessary for the INCUS forward model

This code is comprised of a driver file (postprocess_rams.py) and a file with functions used by the driver (fx_postprocess_rams.py).
The user needs to specify four paths/filenames in the driver file (the rams datafile, the rams headerfile, and the path the save the new data and path to save figures).
The code has been developed for the INCUS mission, and its goal is to provide necessary information for the INCUS forward model.
Specifically, this code calculate 10m wind and surface temperature for the native RAMS variables, and saves these varibles to a new netcdf file.
This code also create a summary figure of the integrated total condensate and vertical velocities in the RAMS data file.
