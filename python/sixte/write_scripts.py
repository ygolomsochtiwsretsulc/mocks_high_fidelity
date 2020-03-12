"""

Writes a number of scripts to be run


"""
import os
import sys
import healpy
import numpy as n

n_pix = healpy.nside2npix(8)

f = open('all_sky_simulate.sh', 'w')
for jj in n.arange(n_pix):
    if jj % 100 == 0:
        f.write('#!/bin/bash \n')
        f.write('source /home/erosita/sw/sass-setup.sh eSASSdevel \n')
    f.write("python simulate.py " + str(jj).zfill(3) + ' \n')

f.close()

#os.system("split -l 102 all_sky_simulate.sh")


f = open('all_sky_pre-process-esass.sh', 'w')
for jj in n.arange(n_pix):
    if jj % 100 == 0:
        f.write('#!/bin/bash \n')
        f.write('source /home/erosita/sw/sass-setup.sh eSASSdevel \n')
    f.write("python pre-process-esass.py " + str(jj).zfill(3) + ' \n')

f.close()


f = open('all_sky_pre-process-esass.sh', 'w')
for jj in n.arange(n_pix):
    if jj % 100 == 0:
        f.write('#!/bin/bash \n')
        f.write('source /home/erosita/sw/sass-setup.sh eSASSdevel \n')
    f.write("python esass.py " + str(jj).zfill(3) + ' \n')

f.close()

for jj in n.arange(n_pix):
    f = open('run_dir/sim_' + str(jj).zfill(3) + '.sh', 'w')
    f.write('#!/bin/bash \n')
    f.write('source /home/erosita/sw/sass-setup.sh eSASSdevel \n')
    f.write("python ../simulate.py " + str(jj).zfill(3) + ' \n')
    #f.write("python ../pre-process-esass.py " + str(jj).zfill(3) + ' \n')
    #f.write("python ../esass.py " + str(jj).zfill(3) + ' \n')
    f.close()

f = open('all_sky_nohup_commands.sh', 'w')
f.write('#!/bin/bash \n')
for jj in n.arange(n_pix):
    f.write(
        'nohup sh sim_' +
        str(jj).zfill(3) +
        '.sh > sim_' +
        str(jj).zfill(3) +
        '.log & \n')

f.close()
