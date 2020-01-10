"""
Creates a fits catalog containing the 4FS input columns.

Create the ID array
remove unwanted column
check rulesets and templates names are correct

"""
import glob
import sys
from astropy_healpix import healpy
import os
from scipy.stats import scoreatpercentile
import pandas as pd  # external package
from scipy.special import erf
from astropy.coordinates import SkyCoord
import astropy.constants as cc
import astropy.io.fits as fits
from astropy.table import Table, Column
import astropy.units as u
import numpy as n
from scipy.stats import norm
from astropy.cosmology import FlatLambdaCDM

cosmoMD = FlatLambdaCDM(
	H0=67.77 * u.km / u.s / u.Mpc,
	Om0=0.307115)  # , Ob0=0.048206)
h = 0.6777
L_box = 1000.0 / h
cosmo = cosmoMD

root_dir = '/home/comparat/wwwDir/eROSITA_AGN_mock/catalogs/MDPL2/'
hd = fits.open(os.path.join(root_dir, 'eRosita_eRASS8_with_photometry.fits'))[1].data
out_dir = '/home/comparat/wwwDir/eROSITA_AGN_mock/catalogs/MDPL2/SC'
quasar = (hd['AGN_type'] == 11)
bright = (hd['HSC-r'] < 21.5) 
s1 = (bright)&(quasar)
N_agn = len(s1.nonzero()[0])
LX = 10**hd['LX_soft'][s1] # erg/s
LX_bin = 2-0.5 # keV
keV_2_Hz = 4.135667662 * 10**(-18) # keV s
#u.keV.to(u.Hz)
LX_monochromatic = n.log10(LX * keV_2_Hz / LX_bin )# erg s-1 Hz-1

def mak_cat(gamma_str = '0.6', scatter_str = '0.2'):
	beta = 8.
	scatter = float(scatter_str)
	gamma = float(gamma_str)
	out_name = os.path.join(out_dir, 'eRosita_eRASS8_with_photometry_G_'+gamma_str+'_S_'+scatter_str+'.fits')
	#LUV_monochromatic_mean = (LX_monochromatic - beta ) / gamma
	#sigma = norm.rvs(loc=0, scale=scatter, size=N_agn)
	#LUV_monochromatic = sigma +  LUV_monochromatic_mean
	sigma = norm.rvs(loc=0, scale=scatter, size=N_agn)
	LUV_monochromatic = (LX_monochromatic - beta - sigma) / gamma
	dL = cosmo.luminosity_distance(hd['redshift_S'][s1]).to(u.cm)
	sphere = (1+hd['redshift_S'][s1]) /(4. * n.pi* (dL.value)**2.)
	FX_monochromatic = 10**LX_monochromatic * sphere
	FUV_monochromatic = 10**LUV_monochromatic * sphere
	t = Table()
	t['RA'] = Column(hd['RA'][s1], unit='degree', dtype=n.float64)
	t['DEC'] = Column(hd['DEC'][s1], unit='degree', dtype=n.float64)
	t['Z'] = Column(hd['redshift_S'][s1], unit='', dtype=n.float32)
	t['FX'] = Column(FX_monochromatic, unit='erg/cm2/s/Hz', dtype=n.float32)
	t['FUV'] = Column(FUV_monochromatic, unit='erg/cm2/s/Hz', dtype=n.float32)
	t['r_mag'] = Column(hd['HSC-r'][s1], unit='mag', dtype=n.float32)
	if os.path.isfile(out_name):
		os.system("rm " + out_name)
	t.write(out_name, format='fits')

mak_cat('0.5' , '0.1')
mak_cat('0.55', '0.1')
mak_cat('0.6' , '0.1')
mak_cat('0.65', '0.1')
mak_cat('0.7' , '0.1')
mak_cat('0.75', '0.1')

mak_cat('0.5' , '0.15')
mak_cat('0.55', '0.15')
mak_cat('0.6' , '0.15')
mak_cat('0.65', '0.15')
mak_cat('0.7' , '0.15')
mak_cat('0.75', '0.15')

mak_cat('0.5' , '0.2')
mak_cat('0.55', '0.2')
mak_cat('0.6' , '0.2')
mak_cat('0.65', '0.2')
mak_cat('0.7' , '0.2')
mak_cat('0.75', '0.2')

mak_cat('0.5' , '0.25')
mak_cat('0.55', '0.25')
mak_cat('0.6' , '0.25')
mak_cat('0.65', '0.25')
mak_cat('0.7' , '0.25')
mak_cat('0.75', '0.25')

mak_cat('0.5' , '0.3')
mak_cat('0.55', '0.3')
mak_cat('0.6' , '0.3')
mak_cat('0.65', '0.3')
mak_cat('0.7' , '0.3')
mak_cat('0.75', '0.3')

mak_cat('0.5' , '0.4')
mak_cat('0.55', '0.4')
mak_cat('0.6' , '0.4')
mak_cat('0.65', '0.4')
mak_cat('0.7' , '0.4')
mak_cat('0.75', '0.4')

mak_cat('0.5' , '0.5')
mak_cat('0.55', '0.5')
mak_cat('0.6' , '0.5')
mak_cat('0.65', '0.5')
mak_cat('0.7' , '0.5')
mak_cat('0.75', '0.5')

