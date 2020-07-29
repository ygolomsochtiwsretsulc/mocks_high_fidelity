"""
What it does
------------

python2.7 on ds43 !!

- 1. Create the cluster X-ray spectra

References
----------

Command to run
--------------

python3 004_7_cluster_cluster_spectra.py

arguments
---------


Dependencies
------------

import xspec

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

"""
import xspec
import sys
import os
import time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import numpy as n
print('CREATES SIMPUT CLUSTER SPECTRA')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits

env = sys.argv[1]  # 'MD04'

# simulation setup
if env[:2] == "MD" : # env == "MD04" or env == "MD40" or env == "MD10" or env == "MD25"
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env[:4] == "UNIT" : # == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

root_dir = os.path.join(os.environ[env])

dir_2_SMPT = os.path.join(root_dir, "cat_CLU_SIMPUT")
dir_2_eRO_all = os.path.join(dir_2_SMPT, "cluster_Xspectra")

if os.path.isdir(dir_2_SMPT) == False:
    os.system('mkdir -p ' + dir_2_SMPT)
if os.path.isdir(dir_2_eRO_all) == False:
    os.system('mkdir -p ' + dir_2_eRO_all)


# NOW links to the grid of SPECTRA
# 0.01 bins for athena
# 0.05 bins for athena
# 0.2 bins for erosita
# kt_arr = n.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
# z_arr = n.hstack((n.array([0., 0.05]), n.arange(0.1, 4., 0.1)))
kt_arr = 10**n.arange(-1,1.3,0.01)
z_arr = n.hstack((n.array([0., 0.01]), 10**n.arange(n.log10(0.02), n.log10(4.), 0.01)))

kt_arrays, z_arrays = n.meshgrid(kt_arr, z_arr)

kt_array = n.hstack((kt_arrays))
z_array = n.hstack((z_arrays))

xspec.Xset.cosmo = str(cosmo.H0.value) + " " + str(cosmo.Ode0)

metal_abundance = 0.3
d_kev = 0.01
#kevs = 10**n.arange(-1., n.log10(50), d_kev)
norm1 = 1.


def get_spectrum(
        temperature=1,
        redshift=0.5,
        metal_abundance=0.3,
        norm1=1.,
        d_kev=0.0025):
    m1 = xspec.Model("apec")
    m1.setPars(
        temperature,  # 1    1   apec       kT         keV      1.00000      +/-  0.0
        metal_abundance,  # 2    1   apec       Abundanc            1.00000      frozen
        redshift,  # 3    1   apec       Redshift            0.0          frozen
        norm1            # 4    1   apec       norm                1.00000      +/-  0.0
    )
    kevs = 10**n.arange(-1., n.log10(50), d_kev)
    fluxes = []
    nPh = []
    for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1], kevs[1:]):
        xspec.AllModels.calcFlux(
            str(kev_min_erosita_RF) +
            " " +
            str(kev_max_erosita_RF))
        fluxes.append(m1.flux[0])
        nPh.append(m1.flux[3])
    return (kevs[:-1] + kevs[1:]) * \
        0.5, n.array(fluxes), n.array(nPh), -kevs[:-1] + kevs[1:]


def create_fits_file_spectrum(outfile, energies, flux_density):
    n_values = str(len(energies))
    print( energies.shape, flux_density.shape)
    print( n.array([energies], dtype='float'))
    col1 = fits.Column(
        name='ENERGY',
        unit='keV',
        format='1PE(' + n_values + ')',
        array=[energies])
    col2 = fits.Column(
        name='FLUXDENSITY',
        unit='photon/s/cm**2/keV',
        format='1PE(' + n_values + '%MCEPASTEBIN%)',
        array=[flux_density])
    cols = fits.ColDefs([col1, col2])
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.name = 'SPECTRUM'
    hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
    hdu.header['HDUCLAS1'] = 'SPECTRUM'
    hdu.header['HDUVERS'] = '1.1.0'
    hdu.header['RADESYS'] = 'FK5'
    hdu.header['EQUINOX'] = 2000.0
    if os.path.isfile(outfile):
        os.system("rm " + outfile)
    hdu.writeto(outfile, clobber=True)


for shell_id, (temperature, redshift) in enumerate(zip(kt_array, z_array)):
    energies, y, nP, dE = get_spectrum(temperature, redshift, metal_abundance=metal_abundance, norm1=norm1, d_kev=d_kev)
    average_flux = y / dE
    # writes
    out_name = os.path.join(dir_2_eRO_all, 'cluster_spectrum_10000kT_' + str(int(temperature * 10000)).zfill(7) + '_10000z_' + str(int(redshift * 10000)).zfill(7) + '.fits')
    print( out_name )
    create_fits_file_spectrum(out_name, energies, average_flux)
