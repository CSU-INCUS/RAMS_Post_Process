#!/bin/bash

# code to run the plotting script
# this is needed because you can't run multithreaded python correctly
# directly from mpiexec in the pleiades script.
cd /home5/pmarines/INCUS/Plot_Code/git/RAMS_Post_Process 

# we have to activate the conda environment, this code will do that.
eval "$(conda shell.bash hook)"
conda activate plotenv

# run the plotting script
# this should NOT be executed in the background (with &). 
# if you do that, the script will terminate immediately. 
# Point to location of jug
/home5/pmarines/miniconda3/envs/plotenv/bin/jug execute plot_w_itc_II.py
