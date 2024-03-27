#!/bin/csh -f 
@ rank = $1 + $MPT_MPI_RANK
# syntax is jug execute check_files_corrupt_jug.py folder_to_look output_loc
/nobackup/jbukowsk/mambaforge/envs/RAMS_QC/bin/jug execute check_files_corrupt_jug.py /nobackup/jbukowsk/PROD/V1/ARG1.1-R-V1/G1/out/ /nobackup/jbukowsk/PROD/RAMS_QC_Check/output_files/ARG1.1_corrupt_files.txt >> /nobackup/jbukowsk/PROD/RAMS_QC_Check/output_files/run_jug_ARG1.1.txt
