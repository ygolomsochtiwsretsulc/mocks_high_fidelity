"""
export GIT_AGN_MOCK='/home/comparat/software/lss_mock_dev/'
export GIT_XPAINT='/home/comparat/software/XSB_Painting'

# run dir :


"""

import time
import sys
import os
from os.path import join
import glob
import numpy as n

lss_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python')
run_dir = os.path.join(lss_git_dir, 'run_dir')


def writeScript(env, baseName):
    f = open(os.path.join(run_dir, env + '_' + baseName + ".sh"), 'w')
    f.write("#!/bin/bash \n")
    f.write("#SBATCH --time=2000:00:00 \n")
    f.write("#SBATCH --nodes=1 \n")
    f.write("#SBATCH --ntasks=1 \n")
    f.write("#SBATCH --cpus-per-task=8 \n")
    f.write("#SBATCH --job-name=" + env + '_' + baseName + " \n")
    f.write(" \n")
    f.write(". /home_local/4FSOpsim/py36he2srv/bin/activate \n")
    f.write("export OMP_NUM_THREADS=1 \n")
    f.write(" \n")
    f.write("cd " + lss_git_dir + " \n")
    f.write(" \n")
    #f.write("ipython3 001_coordinates.py "+env + ' ' + baseName+" \n")
    #f.write("python3 002_0_galaxy.py "+env + ' ' + baseName+" \n")
    #f.write("python3 003_0_agn.py "+env + ' ' + baseName+" \n")
    #f.write("python3 003_1_agn_catalogs.py "+env + ' ' + baseName+" \n")
    f.write(" \n")
    f.close()


t0 = time.time()


env = 'MD04'
test_dir = os.path.join(os.environ[env], 'fits')
baseNames = sorted(n.array([os.path.basename(el)[
                   :-5] for el in n.array(glob.glob(os.path.join(test_dir, '???_?.?????.fits')))]))
for bn in baseNames:
    writeScript(env, bn)


env = 'MD10'
test_dir = os.path.join(os.environ[env], 'hlist', 'fits')
baseNames = n.array([os.path.basename(el)[
                    :-5] for el in n.array(glob.glob(os.path.join(test_dir, '???_?.?????.fits')))])
baseNames.sort()
for bn in baseNames:
    writeScript(env, bn)
