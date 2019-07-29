"""
Author
------

mainly: T. Liu

What it does
------------

python2.7 on ds43 !!

- 1. Create the AGN X-ray spectra

References
----------

Command to run
--------------

python3 003_5_agn_create_spectra.py

arguments
---------


Dependencies
------------

import xspec

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

"""
import sys
import os
import time
import xspec as XS
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import numpy as n
print('CREATES SIMPUT AGN SPECTRA')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

env = sys.argv[1]  # 'MD04'

# simulation setup
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

root_dir = os.path.join(os.environ[env])
dir_2_SMPT = os.path.join(root_dir, "cat_AGN_SIMPUT")
dir_2_eRO_all = os.path.join(dir_2_SMPT, "agn_Xspectra")

if os.path.isdir(dir_2_SMPT) == False:
    os.system('mkdir -p ' + dir_2_SMPT)
if os.path.isdir(dir_2_eRO_all) == False:
    os.system('mkdir -p ' + dir_2_eRO_all)

XS.Xset.cosmo = str(cosmo.H0.value) + " " + str(cosmo.Ode0)


def saveintofile(FileName, N=0):
    if N > 0:
        XS.AllData.dummyrsp(lowE=0.1, highE=50, nBins=N)
    XS.Plot('model')
    Col_E = fits.column.Column(array=n.array(
        [XS.Plot.x(), ], dtype=n.object), name='ENERGY', format='PE()', unit='keV')
    Col_F = fits.column.Column(
        array=n.array(
            [
                XS.Plot.model(),
            ],
            dtype=n.object),
        name='FLUXDENSITY',
        format='PE()',
        unit='photon/s/cm**2/keV')
    HDU1 = fits.BinTableHDU.from_columns([Col_E, Col_F])
    HDU1.header['EXTNAME'] = 'SPECTRUM'
    HDU1.header['HDUCLAS1'] = 'SPECTRUM'
    HDU1.writeto(os.path.join(dir_2_eRO_all, FileName), overwrite=True)


XS.Model('TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)')
pars = ['0.01    -0.0001          0          0     100000      1e+06',
        '1       0.01      1e-06      1e-06        1e5        1e5',
        '3',
        '1      -0.01          0          0         10         10',
        '7.11    -0.0711          7          7         10         10',
        '1.9     -0.019          0          0          3          3',
        '95      -0.95       0.01          1        100        200',
        '300         -3          1          1      1e+06      1e+06',
        '1      -0.01          0          0          1          1',
        '0',
        '0.5     -0.005     -0.999     -0.999         10         10',
        '1      -0.01          0          0      1e+20      1e+24',
        '6.4     -0.064          0          0      1e+06      1e+06',
        '0.05    -0.0005          0          0         10         20',
        '= p11',
        '0.01 2.6316e-11          0          0      1e+20      1e+24',
        '0.02    -0.0002          0          0        0.1        0.1',
        '= p6',
        '= p12/(1. + p11)/(1./(1. + p11))^( - p6)',
        '= p6',
        '300         -3          1          1      1e+06      1e+06',
        '-1      -0.01         -3         -3     -1e-09     -1e-09',
        '= p11',
        '1      -0.01          0          0      1e+06      1e+06',
        '1      -0.01          0          0      1e+06      1e+06',
        '0.45    -0.0045       0.05       0.05       0.95       0.95',
        '= p19',
        '1      -0.01       0.01       0.01        100        100']
XS.AllModels(1).setPars(pars)

for nH in n.arange(20, 26.2, 0.2):
    XS.AllModels(1).plcabs.nH = 10**(nH - 22)
    for z in n.arange(0, 6.1, 0.1):
        XS.AllModels(1).plcabs.Redshift = z
        for nb in 2**n.arange(2, 11):
            filename = 'NH' + str(n.round(nH, 1)) + '_Z' + \
                str(n.round(z, 1)) + '_N' + str(int(nb)) + '.fits'
            print(filename)
            saveintofile(filename, nb)
