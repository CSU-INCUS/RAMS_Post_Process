#Sean.Freeman@colostate.edu
#python script to check for corrupt RAMS files
#usage: jug execute check_files_corrupt_jug.py folder_to_look output_loc
#where type is 'analysis' or 'lite'
#e.g. jug execute check_files_corrupt_jug.py ./ corrupt_file_log.txt

import sys

if len(sys.argv) != 3:
    print("Usage: jug execute check_files_corrupt_jug.py folder_to_look log_output")
    sys.exit()

import h5py
import numpy
import glob
import hdf5plugin
from jug import TaskGenerator

glob_look = sys.argv[1]+'/a-?-*.h5'

@TaskGenerator
def check_file_corrupt(in_file_name):
    test_in = h5py.File(in_file_name, 'r')
    print(in_file_name)
    corrupt_vars = list()
    for var in test_in:
        try:
            x = test_in[var][0]
        except:
            corrupt_vars.append(var)
    return corrupt_vars

@TaskGenerator
def write_corrupt_files(corrupt_file_dict, outloc):
    with open(outloc, 'w') as outfil:
        for entry in corrupt_file_dict:
            if len(corrupt_file_dict[entry])>0:
                print(entry+" corrupt on these variables: "+" ".join(corrupt_file_dict[entry]), file=outfil)
                        
corrupt_file_dict = dict()
for filename in glob.glob(glob_look):
    outfilcor = check_file_corrupt(filename)
    corrupt_file_dict[filename] = outfilcor

write_corrupt_files(corrupt_file_dict, sys.argv[-1])
