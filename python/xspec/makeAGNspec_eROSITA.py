#!/usr/bin/env python2
"""
Authors: T. Liu (MPE/MPG)
Date: 11.03.2019

Generates a set of representative templates for eRosita AGNs.

"""

import os
import sys
import numpy as np
from astropy.io import fits
import xspec as XS


def saveintofile(FileName, N=0):
    if N > 0:
        XS.AllData.dummyrsp(lowE=0.1, highE=50, nBins=N)
    XS.Plot('model')
    Col_E = fits.column.Column(array=np.array(
        [XS.Plot.x(), ], dtype=np.object), name='ENERGY', format='PE()', unit='keV')
    Col_F = fits.column.Column(
        array=np.array(
            [
                XS.Plot.model(),
            ],
            dtype=np.object),
        name='FLUXDENSITY',
        format='PE()',
        unit='photon/s/cm**2/keV')
    HDU1 = fits.BinTableHDU.from_columns([Col_E, Col_F])
    HDU1.header['EXTNAME'] = 'SPECTRUM'
    HDU1.header['HDUCLAS1'] = 'SPECTRUM'
    HDU1.writeto(FileName, overwrite=True)


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

for nH in np.arange(20, 26.2, 0.2):
    print(nH)
    XS.AllModels(1).plcabs.nH = 10**(nH - 22)
    for z in np.arange(0, 6.1, 0.1):
        print(z)
        XS.AllModels(1).plcabs.Redshift = z
        for nb in 2**np.arange(2, 11):
            print(nb)
            tpl_name = '/afs/mpe/www/people/comparat/eROSITA_AGN_mock/spectra/Xray/spectra/' + 'AGN_NH_' + \
                str(np.round(nH, 1)) + '_Z_' + str(np.round(z, 1)) + '_N_' + str(int(nb)) + '.fits'
            print(tpl_name)
            saveintofile(tpl_name, nb)
