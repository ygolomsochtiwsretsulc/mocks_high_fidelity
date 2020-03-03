"""
What it does
------------

Makes a small fits file out of the full N-body fits file shell.

command to run
--------------

python3 Create_small_fits_for_plotting.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

topcat / stilts

"""

import numpy as n
import os
import sys
import glob

env = sys.argv[1]  # 'MD04'

top_dir = os.path.join(os.environ[env] )
fits_dir = os.path.join(top_dir, 'fits')

out_dir = os.path.join(top_dir, 'mini_fits')
if os.path.isdir(out_dir) == False:
    os.mkdir(out_dir)

fits_list = sorted(
    n.array(
        glob.glob(
            os.path.join(
                fits_dir,
                'all_?.?????.fits'))))

for in_file in fits_list[::-1]:
    out_file = os.path.join(out_dir, os.path.basename(in_file))
    stilts_command = "stilts tpipe in=" + in_file + """ ifmt=fits cmd='select "z>-5 && z<5"' omode=out ofmt=fits out=""" + out_file
    print(stilts_command)
    os.system(stilts_command)

fits_list = n.array(glob.glob(os.path.join(fits_dir, 'sat_?.?????.fits')))
fits_list.sort()

for in_file in fits_list[::-1]:
    out_file = os.path.join(out_dir, os.path.basename(in_file))
    stilts_command = "stilts tpipe in=" + in_file + """ ifmt=fits cmd='select "z>-5 && z<5"' omode=out ofmt=fits out=""" + out_file
    print(stilts_command)
    os.system(stilts_command)
