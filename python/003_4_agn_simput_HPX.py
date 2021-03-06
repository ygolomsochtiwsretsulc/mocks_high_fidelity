"""
What it does
------------

Creates a simput file AGN fits catalog for each healpix 768 pixel of 13 deg2 (NSIDE=8)

Command to run
--------------

python3 003_1_agn_catalogs.py environmentVAR f_sat laptop

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

f_sat: fraction of satellite e.g. 0.1

laptop: if running on the server or on the laptop. Redirects properly to the stilts command

Dependencies
------------
topcat/stilts
import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib

"""
from astropy_healpix import healpy
import sys
import os
import time
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import h5py
import numpy as n
print('CREATES FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits
# import all pathes

env = sys.argv[1]  # 'MD04'
HEALPIX_id = int(sys.argv[2])  # 700
print(sys.argv) 

stilts_cmd = 'stilts'

root_dir = os.path.join(os.environ[env])

dir_2_eRO_all = os.path.join(root_dir, "cat_AGN_all")

dir_2_SMPT = os.path.join(root_dir, "cat_AGN_SIMPUT")

if os.path.isdir(dir_2_SMPT) == False:
    os.system('mkdir -p ' + dir_2_SMPT)

path_2_eRO_all_catalog = os.path.join(dir_2_eRO_all, str(HEALPIX_id).zfill(6) + '.fit')
path_2_SMPT_catalog = os.path.join(dir_2_SMPT, 'SIMPUT_' + str(HEALPIX_id).zfill(6) + '.fit')
print('=================================================================')
print(path_2_eRO_all_catalog, path_2_SMPT_catalog) # path_2_eRO_sat_catalog

hd_all = fits.open(path_2_eRO_all_catalog)
N_agn_all = len(hd_all[1].data['ra'])
N_all = int(N_agn_all)# * (1 - f_sat))

indexes_all = (n.arange(N_agn_all) + 1e9 ).astype('int')

rd_all = n.random.rand(N_agn_all)

# FLUX limit in the X-ray
FX_LIM_value_cen = 2e-17 
detected_all = (hd_all[1].data['AGN_FX_soft'] > FX_LIM_value_cen)

# overall selection of central and satellites
sel_all = (detected_all) #  (rd_all < 1 - f_sat) &

# stacking dara arrays 
# indexes
data_indexes = indexes_all[sel_all] 
# indexes_all
# NH
data_nh = hd_all[1].data['AGN_Nh'][sel_all] 
data_nh = (data_nh * 5).astype('int') / 5.
data_nh[data_nh < 20.] = 20.0
data_nh[data_nh > 26.] = 26.0

data_z = n.round(hd_all[1].data['redshift_R'][sel_all],1) 
ra_array = hd_all[1].data['ra'][sel_all] 
dec_array = hd_all[1].data['dec'][sel_all] 

n_e_bins = 2**n.arange(2, 11)

n_total = int(512e6) 
n_allowed = (n_total / n_e_bins).astype('int')[::-1] - 100

FX_array = hd_all[1].data['AGN_FX_soft'][sel_all] 

fbins = n.arange(-n.log10(FX_array).max() - 0.1, -
					n.log10(FX_array).min() + 0.1, 0.01)
xf = fbins[:-1] + 0.01 / 2
hst = n.cumsum(n.histogram(-n.log10(FX_array), bins=fbins)[0])
itp = interp1d(hst, xf)

n_allowed_c = n.cumsum(n_allowed)
n_allowed_t = n_allowed_c[n_allowed_c < itp.x[-1]]
FX_inner_boundaries = -itp(n_allowed_t)
FX_boundaries = n.hstack((
	n.log10(FX_array).max(),
	FX_inner_boundaries,
	n.log10(FX_array).min()))

data_n_e_b = (n.ones(len(data_z)) * n_e_bins[-1]).astype('int')

for jj, (f_min, f_max) in enumerate(
		zip(FX_boundaries[:-1], FX_boundaries[1:])):
	selection = (n.log10(FX_array) <= f_min) & (n.log10(FX_array) >= f_max)
	n_e_val = n_e_bins[::-1][jj]
	data_n_e_b[selection] = n_e_val
	print(f_min,
			f_max,
			n_e_val,
			len(data_n_e_b[selection]) < 512e6 / n_e_val,
			len(data_n_e_b[selection]),
			512e6 / n_e_val,
			data_n_e_b[selection])

# NH24.2_Z3.9_1024.fits
# 'NH'+str(n.round(nH,1))+'_Z'+str(n.round(z,1))+'_N'+str(int(nb))+'.fits'
s1 = 'agn_Xspectra/NH'
str_nh = n.array(data_nh.astype('str'), dtype=n.object)
s2 = '_Z'
str_z = n.array(data_z.astype('str'), dtype=n.object)
s3 = '_N'
str_n_e_b = n.array(data_n_e_b.astype('str'), dtype=n.object)
s4 = '.fits' + """[SPECTRUM][#row==1]"""

tpl = s1 + str_nh + s2 + str_z + s3 + str_n_e_b + s4
print(tpl)
# /afs/mpe/www/people/comparat/eROSITA_AGN_mock/spectra/Xray/spectra/
# 'agn/spectra/agn_nH_21.6_z_4.1_nEbins_539.fits[SPECTRUM][#row==1]'

N_agn_out = len(ra_array)

hdu_cols = fits.ColDefs([
	fits.Column(name="SRC_ID", format='K', unit='', array=data_indexes), 
	fits.Column(name="RA", format='D', unit='deg', array=ra_array), 
	fits.Column(name="DEC", format='D', unit='deg', array=dec_array), 
	fits.Column(name="E_MIN", format='D', unit='keV', array=n.ones(N_agn_out) * 0.5), 
	fits.Column(name="E_MAX", format='D', unit='keV', array=n.ones(N_agn_out) * 2.0), 
	fits.Column(name="FLUX", format='D', unit='erg/s/cm**2', array=FX_array), 
	fits.Column(name="SPECTRUM", format='100A', unit='', array=tpl), 
	fits.Column(name="n_energy_bins", format='K', unit='', array=data_n_e_b)
])

hdu = fits.BinTableHDU.from_columns(hdu_cols)

hdu.name = 'SRC_CAT'
hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
hdu.header['HDUCLAS1'] = 'SRC_CAT'
hdu.header['HDUVERS'] = '1.1.0'
hdu.header['RADESYS'] = 'FK5'
hdu.header['EQUINOX'] = 2000.0

outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])
if os.path.isfile(path_2_SMPT_catalog):
	os.system("rm " + path_2_SMPT_catalog)
outf.writeto(path_2_SMPT_catalog, overwrite=True)
print(path_2_SMPT_catalog, 'written', time.time() - t0)

for neb in n_e_bins[::-1][:1]:
	path_2_SMPT_catalog_slice = os.path.join(
		dir_2_SMPT,
		'SIMPUT_' +
		str(HEALPIX_id).zfill(6) +
		'_{}.fit'.format(neb))
	print('slice after n e bins', time.time() - t0)
	#c1="java -jar ~/software/stilts.jar tpipe ifmt=fits in="+path_2_SMPT_catalog
	c1 = stilts_cmd + " tpipe ifmt=fits in=" + path_2_SMPT_catalog
	c2 = """ cmd='select "n_energy_bins=={}"' """.format(neb)
	c3 = " omode=out ofmt=fits out={}".format(path_2_SMPT_catalog_slice)
	print(c1 + c2 + c3)
	os.system(c1 + c2 + c3)

os.system('rm ' + path_2_SMPT_catalog)

#print('check for missing templates')
# check that all spectra are there
#hh = fits.open(path_2_SMPT_catalog)
#all_specs = n.unique(hh[1].data['SPECTRUM'])
#spec_dir = '/afs/mpe/www/people/comparat/eROSITA_AGN_mock/spectra/Xray/'
# for spec in all_specs:
# if os.path.isfile(os.path.join(spec_dir, spec[:-19]))==False:
#print(os.path.join(spec_dir, spec[:-19]), spec)
