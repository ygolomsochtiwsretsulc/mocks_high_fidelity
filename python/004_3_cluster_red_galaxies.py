"""
What it does
------------

Take the complete cluster galaxy satellite FILE:  $ENV_eRO_CLU_SAT.fit
Re-assigns the is_quiescent flag as a function of radius as well as star formation rates.

It increases the number of quiescent galaxies in the vicinity of clusters

References
----------

 * Hennig et al. 2017, https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.4015H
 * Biffi et al. Dolag et al. Magneticum, https://ui.adsabs.harvard.edu/#abs/2018MNRAS.tmp.2317B

Command to run
--------------

python3 004_3_cluster_red_galaxies.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, numpy, astropy


"""

from scipy.optimize import newton
import sys
import os
import time
from scipy.stats import norm
import astropy.io.fits as fits
import numpy as n
print('Adjusts red sequence of galaxies around clusters')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

# import all pathes
env = sys.argv[1]

test_dir = os.path.join(os.environ[env], 'hlists', 'fits')
path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT.fit')


# INCREASE QUIESCENT FRACTIONS !

hdu_clu = fits.open(path_2_CLU_SAT_catalog, mode='update')
x = hdu_clu[1].data['comoving_distance_to_cluster_in_rvir']
is_quiescent = hdu_clu[1].data['is_quiescent']
zr_CLU = hdu_clu[1].data['redshift_R']
mass = hdu_clu[1].data['galaxy_stellar_mass']
log_sfr = hdu_clu[1].data['galaxy_star_formation_rate']


def frac_old(x, z_cluster): return ((x + 0.01)**(-0.25) - \
             x / 100. - 0.47) * (1. + z_cluster)**2. / 3.2


f_red = frac_old(x, zr_CLU)

rds = n.random.random(len(x))

is_red = (rds < f_red)
is_1rvir = (x < 1)

print(len(is_red[is_1rvir].nonzero()[0]))
print(len(is_quiescent[is_1rvir].nonzero()[0]))

N_red_1rvir_needed = len(is_red[is_1rvir].nonzero()[0])
N_quiescent_1rvir_already_there = len(is_quiescent[is_1rvir].nonzero()[0])

N_red_1rvir_to_assign = N_red_1rvir_needed - N_quiescent_1rvir_already_there

# find them among the non_quiescent
SF = (is_quiescent == False)
rds2 = n.random.random(len(x[SF]))
is_1rvir2 = (x[SF] < 1)
f_red2 = frac_old(x[SF], zr_CLU[SF])


def func(DELTA_X):
    change_2_red = (rds2 < f_red2 - DELTA_X)
    N_new_red = len(change_2_red[is_1rvir2].nonzero()[0])
    # print(N_new_red)
    return N_red_1rvir_to_assign - N_new_red  # , change_2_red


VAL = newton(func, 0.1)
change_2_red = (rds2 < f_red2 - VAL)[(x[SF] < 1)]

is_quiescent[(SF) & (x < 1)] = change_2_red

# change the SFR values
# mass-SFR sequence for the quiescent
sfr_Q = n.zeros_like(zr_CLU)


def beta_z(z): return -0.57 * z + 1.43


def alpha_z(z): return 6.32 * z - 16.26


def mean_SFR_Q(mass, z): return mass * beta_z(z) + alpha_z(z)


def scale_z(z): return -0.34 * z + 0.99


rds3 = norm.rvs(loc=0, scale=1., size=len(
    zr_CLU[(SF) & (x < 1)])) * scale_z(zr_CLU[(SF) & (x < 1)])
log_sfr_Q = mean_SFR_Q(mass[(SF) & (x < 1)], zr_CLU[(SF) & (x < 1)]) + rds3
# change SFR for the quiesent selection
log_sfr[(SF) & (x < 1)] = log_sfr_Q


# update the file SFR and is_quiescent columns
hdu_clu[1].data['is_quiescent'] = is_quiescent
hdu_clu[1].data['galaxy_star_formation_rate'] = log_sfr

hdu_clu.flush()
