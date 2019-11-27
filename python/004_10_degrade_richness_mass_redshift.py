"""
What it does
------------

Creates a mock with degrade redshifts, masses and richnesses 

"""

import glob
import sys
from astropy_healpix import healpy
import os
import time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import healpy as hp
from scipy.stats import norm
from scipy.interpolate import interp1d
import astropy.io.fits as fits
from astropy.table import Table
import numpy as n
print('Degrades a mock catalogue to input cosmopipes')
print('------------------------------------------------')
t0 = time.time()

# simulation name
env = 'MD10' # sys.argv[1]

# cosmology set up
if env == "MD10" or env == "MD04":
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

# where the catalogue is :
test_dir = os.path.join( os.environ[env] )
# catalog
path_2_RS_catalog = os.path.join(test_dir, 'cat_eRO_CLU_RS', '000356.fit')

tt = Table.read(path_2_RS_catalog)

# degrade the catalogue here
tt['richness'] = tt['richness']/2.

# save the new catalogue here :
path_2_out_file = os.path.join(test_dir, env+'degraded_mock_Lambda_div_2.fits')

tt.write(path_2_out_file, format='fits')



