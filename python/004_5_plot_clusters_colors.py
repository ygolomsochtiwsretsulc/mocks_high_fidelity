"""
What it does
------------

Plots the cluster model

References
----------

Command to run
--------------

python3 004_5_plot_clusters.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then	 work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

"""
from astropy.table import Table, Column
from scipy.stats import scoreatpercentile
import glob
import sys
from astropy_healpix import healpy
import os
import time
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import numpy as n
print('Create file with galaxies around clusters')
print('=> Abundance matching for magnitudes')
print('=> Red sequence colors')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


env =  sys.argv[1]
#delta_crit = sys.argv[2]
#print(env, delta_crit)
test_dir = os.path.join(os.environ[env])
all_catalogs = n.array(glob.glob(os.path.join(test_dir, 'fits', 'all_*_galaxiesAroundClusters.fit')))
all_aexp = n.array([ os.path.basename(el).split('_')[1] for el in all_catalogs ])
all_aexp.sort()

fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'clusters', 'galaxies')
if os.path.isdir(fig_dir) == False:
	os.system('mkdir -p ' + fig_dir)

path_2_spiders_host_file = os.path.join(
	os.environ['HOME'],
	'hegcl/SPIDERS',
	'mastercatalogue_FINAL_CODEXID.fits')
spiders_host = fits.open(path_2_spiders_host_file)[1].data

path_2_spiders_file = os.path.join(
	os.environ['HOME'],
	'hegcl/SPIDERS',
	'mastercatalogue_FINAL_CODEXID-flat_DR16_firefly.fits')
spiders_gal = fits.open(path_2_spiders_file)[1].data

magi_low, magi_high, frac_obs = n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], 'data/spiders/fraction-observed-cmodemag-i.txt'), unpack=True)
frac_obs_itp = interp1d (
	n.hstack(( magi_low[0]-10,(magi_low+ magi_high)/2., magi_high[-1]+10 )) , 
	n.hstack((1., frac_obs, 1. )) ) 


# red sequence DATA
red_seq_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'redSequence')
path_2_f1 = os.path.join(red_seq_dir, 'sdss_cal_iter1_pars.fit')
path_2_f2 = os.path.join(red_seq_dir, 'sdss_cal_zspec_redgals_model.fit')

f1 = fits.open(path_2_f1)
f2 = fits.open(path_2_f2)

redshift_RS = f2[1].data['nodes'][0]
ug_RS = f2[1].data['meancol'].T[0].T[0]
gr_RS = f2[1].data['meancol'].T[1].T[0]
ri_RS = f2[1].data['meancol'].T[2].T[0]
iz_RS = f2[1].data['meancol'].T[3].T[0]

ug_sigma_RS = f2[1].data['meancol_scatter'].T[0].T[0]
gr_sigma_RS = f2[1].data['meancol_scatter'].T[1].T[0]
ri_sigma_RS = f2[1].data['meancol_scatter'].T[2].T[0]
iz_sigma_RS = f2[1].data['meancol_scatter'].T[3].T[0]

ug_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
					n.hstack((ug_RS[0], ug_RS, ug_RS[-1])))
gr_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
					n.hstack((gr_RS[0], gr_RS, gr_RS[-1])))
ri_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
					n.hstack((ri_RS[0], ri_RS, ri_RS[-1])))
iz_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
					n.hstack((iz_RS[0], iz_RS, iz_RS[-1])))

ug_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
	(ug_sigma_RS[0], ug_sigma_RS, ug_sigma_RS[-1])))
gr_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
	(gr_sigma_RS[0], gr_sigma_RS, gr_sigma_RS[-1])))
ri_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
	(ri_sigma_RS[0], ri_sigma_RS, ri_sigma_RS[-1])))
iz_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
	(iz_sigma_RS[0], iz_sigma_RS, iz_sigma_RS[-1])))


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

# opti0nal red sequence figure
fig_out = os.path.join(fig_dir, 'red_sequence.png')
color_bins = n.arange(-0.0,1.2,0.01)
z_all = n.arange(0., 1.2, 0.05)
p.figure(1, (6., 6.))

p.plot(z_all, gr_RS_itp(z_all), label='gr model', color='blue', ls='solid')
p.plot(z_all, gr_RS_itp(z_all)+gr_sigma_RS_itp(z_all), color='blue', ls='dashed')
p.plot(z_all, gr_RS_itp(z_all)-gr_sigma_RS_itp(z_all), color='blue', ls='dashed')

p.plot(z_all, ri_RS_itp(z_all), label='ri model', color='green', ls='solid')
p.plot(z_all, ri_RS_itp(z_all)+ri_sigma_RS_itp(z_all), color='green', ls='dashed')
p.plot(z_all, ri_RS_itp(z_all)-ri_sigma_RS_itp(z_all), color='green', ls='dashed')

p.plot(z_all, iz_RS_itp(z_all), label='iz model', color='red', ls='solid')
p.plot(z_all, iz_RS_itp(z_all)+iz_sigma_RS_itp(z_all), color='red', ls='dashed')
p.plot(z_all, iz_RS_itp(z_all)-iz_sigma_RS_itp(z_all), color='red', ls='dashed')

p.legend(frameon=False)
p.xlabel('redshift')
p.ylabel('color')
#p.ylabel('probability distribution function')
p.grid()
# p.ylim((0,1.1))
# p.yscale('log')
# p.xscale('log')
p.savefig(fig_out)
p.clf()


for aexp_str in all_aexp[::-1]:
	baseName = 'all_'+aexp_str # sys.argv[2]  # "all_0.62840"
	z_snap = 1./float(baseName.split('_')[1])-1.
	aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
	print(env, baseName)
	test_dir = os.path.join(os.environ[env])

	# import all pathes
	path_2_CLU_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')
	path_2_CLU_SAT_catalog = os.path.join(test_dir, 'fits', baseName + '_galaxiesAroundClusters.fit')
	#path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT.fit')
	#path_2_CLU_SAT_RS_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT_RS.fit')

	hdu_clu = Table.read(path_2_CLU_SAT_catalog)  # , mode='update')
	#x = ((hdu_clu['x']-hdu_clu['HOST_HALO_x'])**2. + (hdu_clu['y']-hdu_clu['HOST_HALO_y'])**2. + (hdu_clu['z']-hdu_clu['HOST_HALO_z'])**2.)**0.5
	is_quiescent = hdu_clu['is_quiescent']
	zr_gal = hdu_clu['redshift_R']
	mean_z = n.mean(zr_gal)
	mag_available = (hdu_clu['sdss_r_err'] > 0 ) & ( is_quiescent )
	gr = hdu_clu['sdss_g'][mag_available] - hdu_clu['sdss_r'][mag_available]
	ri = hdu_clu['sdss_r'][mag_available] - hdu_clu['sdss_i'][mag_available]
	iz = hdu_clu['sdss_i'][mag_available] - hdu_clu['sdss_z'][mag_available]
	
	# opti0nal red sequence figure
	fig_out = os.path.join(fig_dir, 'red_sequence_'+aexp_str+'.png')
	color_bins = n.arange(-0.0,1.2,0.01)
	p.figure(1, (6., 6.))
	p.hist(gr, bins=color_bins, histtype = 'step', rasterized = True,label='gr mock', color='blue' , cumulative=True, density=True, lw=3)
	p.hist(ri, bins=color_bins, histtype = 'step', rasterized = True,label='ri mock', color='green', cumulative=True, density=True, lw=3)
	p.hist(iz, bins=color_bins, histtype = 'step', rasterized = True,label='iz mock', color='red'  , cumulative=True, density=True, lw=3)
	p.axvline(gr_RS_itp(mean_z), label='gr model', color='blue', ls='dotted')
	p.axvline(gr_RS_itp(mean_z) + gr_sigma_RS_itp(mean_z), color='blue', ls='dashed')
	p.axvline(gr_RS_itp(mean_z) - gr_sigma_RS_itp(mean_z), color='blue', ls='dashed')

	p.axvline(ri_RS_itp(mean_z), label='ri model', color='green', ls='dotted')
	p.axvline(ri_RS_itp(mean_z) + ri_sigma_RS_itp(mean_z), color='green', ls='dashed')
	p.axvline(ri_RS_itp(mean_z) - ri_sigma_RS_itp(mean_z), color='green', ls='dashed')

	p.axvline(iz_RS_itp(mean_z), label='iz model', color='red', ls='dotted')
	p.axvline(iz_RS_itp(mean_z) + iz_sigma_RS_itp(mean_z), color='red', ls='dashed')
	p.axvline(iz_RS_itp(mean_z) - iz_sigma_RS_itp(mean_z), color='red', ls='dashed')
	
	#p.errorbar(redshift_RS, gr_RS, xerr=0.05, yerr=gr_sigma_RS, label='g-r')
	#p.errorbar(redshift_RS, ri_RS, xerr=0.05, yerr=ri_sigma_RS, label='r-i')
	#p.errorbar(redshift_RS, iz_RS, xerr=0.05, yerr=iz_sigma_RS, label='i-z')
	p.legend(frameon=False)
	# p.title('z='+str(z_cluster))
	p.xlabel('color')
	p.ylabel('probability distribution function')
	p.grid()
	# p.ylim((0,1.1))
	# p.yscale('log')
	# p.xscale('log')
	p.savefig(fig_out)
	p.clf()
	
