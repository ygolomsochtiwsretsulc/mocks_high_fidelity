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
print('Plot luminosity functino with SPIDERS data')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

env =  sys.argv[1]
delta_crit = sys.argv[2]
print(env, delta_crit)
test_dir = os.path.join(os.environ[env])
all_catalogs = n.array(glob.glob(os.path.join(test_dir, 'fits', 'all_*_galaxiesAroundClusters.fit')))
all_aexp = n.array([ os.path.basename(el).split('_')[1] for el in all_catalogs ])
all_aexp.sort()

#N_rvir_3D = 2.
#delta_crit = 'all'
#delta_crit = '200c'

# Possibly implement the double schechter from Popesso 2005 ?
logPhi_z, alpha_z, M_s_z = n.loadtxt(os.path.join(
    os.environ['GIT_AGN_MOCK'], 'data', 'cluster_galaxies', 'ricci_2018_table3.ascii'))


# schechter function, Eq. 8,9.
def Schechter_L(L, phi_s, L_s, alpha):
	"""
	Eq. 8 Ricci 2018
	"""
	return phi_s * (L / L_s)**alpha * n.e**(-L / L_s) / L_s

def Schechter_M(M, phi_s, M_s, alpha): 
	"""
	Eq. 9 Ricci 2018
	"""
	return 0.4 * n.log(10.) * phi_s * 10**(0.4 * (M_s - M) * (alpha + 1)) * n.e**(-10**(0.4 * (M_s - M)))

# redshift evolution function for the parameters o fthe Schecter function, Eq. 10
def function_eq10(redshift, richness, a, b, c): 
	"""
	Eq. 10 for a median redshift and a median richness
	"""
	return a * n.log10(1 + redshift) + b * n.log10(richness) + c

def logPhi_evol( redshift, richness):
	"""
	Redshift evolution of Phi star
	Combination of the values in Table 3 and Eq. 10 of Ricci 2018
	"""
	return function_eq10( redshift, richness, logPhi_z[0], logPhi_z[1], logPhi_z[2] )

def alpha_evol( redshift, richness): 
	"""
	Redshift evolution of alpha
	Combination of the values in Table 3 and Eq. 10 of Ricci 2018
	"""
	return function_eq10( redshift, richness, alpha_z[0], alpha_z[1], alpha_z[2] )

def M_s_evol(redshift, richness): 
	"""
	Redshift evolution of M star
	Combination of the values in Table 3 and Eq. 10 of Ricci 2018
	"""
	return function_eq10( redshift, richness, M_s_z[0], M_s_z[1], M_s_z[2] )


def Schechter_M_z(M, redshift, richness): 
	"""
	Schechter function as a function of absolute magnitude for a median redshift and a median richness
	Combination of the values in Table 3 and Eq. 8, 9, 10 of Ricci 2018
	"""
	return 0.4 * n.log(10.) * 10**logPhi_evol(redshift, richness) * 10**(0.4 * (M_s_evol(redshift, richness) - M) * (alpha_evol(redshift, richness) + 1)) * n.e**( -10** ( 0.4 * (M_s_evol(redshift,richness) - M)))


def mass_2_richness(M200c, redshift): 
	"""
	M200c to richness conversion using the scaling relation 
	Table 2, second line of Capasso et al. 2019 (1812.0609
	"""
	return 39.8 * (M200c/3e14)**(0.98) * ((1+redshift)/(1+0.18))**(-1.08)

def Schechter_M_z_M200c(M, redshift, M200c): 
	"""
	Schechter function as a function of absolute magnitude for a median redshift and a median M200c
	Combination of the values in Table 3 and Eq. 8, 9, 10 of Ricci 2018
	"""
	return 0.4 * n.log(10.) * 10**logPhi_evol(redshift, mass_2_richness(M200c, redshift)) * 10**(0.4 * (M_s_evol(redshift, mass_2_richness(M200c, redshift)) - M) * (alpha_evol(redshift, mass_2_richness(M200c, redshift)) + 1)) * n.e**( -10** ( 0.4 * (M_s_evol(redshift,mass_2_richness(M200c, redshift)) - M)))


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

mag_R18_015, phi_up_R18_015, phi_low_R18_015 = n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'Ricci_2018', 'LF_fig12_0_z_025.txt'), unpack=True)
mag_R18_030, phi_up_R18_030, phi_low_R18_030 = n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'Ricci_2018', 'LF_fig12_025_z_035.txt'), unpack=True)

for aexp_str in all_aexp[::-1]:
	baseName = 'all_'+aexp_str # sys.argv[2]  # "all_0.62840"
	z_snap = 1./float(baseName.split('_')[1])-1.
	aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
	print(env, baseName)
	test_dir = os.path.join(os.environ[env])

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

	# import all pathes
	path_2_CLU_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')
	path_2_CLU_SAT_catalog = os.path.join(test_dir, 'fits', baseName + '_galaxiesAroundClusters.fit')
	#path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT.fit')
	#path_2_CLU_SAT_RS_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT_RS.fit')

	"""
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

	# opti0nal red sequence figure

	fig_out = os.path.join(red_seq_dir, 'red_sequence.png')

	p.figure(1, (6., 6.))

	p.errorbar(redshift_RS, gr_RS, xerr=0.05, yerr=gr_sigma_RS, label='g-r')
	p.errorbar(redshift_RS, ri_RS, xerr=0.05, yerr=ri_sigma_RS, label='r-i')
	p.errorbar(redshift_RS, iz_RS, xerr=0.05, yerr=iz_sigma_RS, label='i-z')

	p.legend(frameon=False)
	# p.title('z='+str(z_cluster))
	p.xlabel('redshift')
	p.ylabel('color')
	p.grid()
	# p.ylim((0,1.1))
	# p.yscale('log')
	# p.xscale('log')
	p.savefig(fig_out)
	p.clf()
	"""

	# opens mock catalogs

	# galaxies
	#hdu_clu = fits.open(path_2_CLU_SAT_catalog)  # , mode='update')
	hdu_clu = Table.read(path_2_CLU_SAT_catalog)  # , mode='update')
	x = ((hdu_clu['x']-hdu_clu['HOST_HALO_x'])**2. + (hdu_clu['y']-hdu_clu['HOST_HALO_y'])**2. + (hdu_clu['z']-hdu_clu['HOST_HALO_z'])**2.)**0.5
	is_quiescent = hdu_clu['is_quiescent']
	zr_gal = hdu_clu['redshift_R']
	z_cluster = hdu_clu['HOST_redshift_R']
	DM = hdu_clu['mag_abs_r'] - hdu_clu['mag_r'] 
	mock_mag_abs = hdu_clu['mag_abs_r'] 
	#mock_mag_abs_sdss_u = hdu_clu['sdss_u'] + DM
	#mock_mag_abs_sdss_g = hdu_clu['sdss_g'] + DM
	#mock_mag_abs_sdss_r = hdu_clu['sdss_r'] + DM
	#mock_mag_abs_sdss_i = hdu_clu['sdss_i'] + DM
	#mock_mag_abs_sdss_z = hdu_clu['sdss_z'] + DM

	omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
	DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
	HOST_HALO_Mvir = hdu_clu['HOST_HALO_Mvir'] / h
	HOST_HALO_Rvir = hdu_clu['HOST_HALO_Rvir']
	HOST_HALO_M500c = hdu_clu['HOST_HALO_M500c'] / h
	HOST_HALO_R500c = (DeltaVir_bn98(z_snap)/500. * HOST_HALO_M500c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir
	HOST_HALO_M200c = hdu_clu['HOST_HALO_M200c'] / h
	HOST_HALO_R200c = (DeltaVir_bn98(z_snap)/200. * HOST_HALO_M200c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir

	zmin = n.min(zr_gal) 
	zmax = n.max(zr_gal) 
	print(zmin, '<z<', zmax)
	if delta_crit == '200c' :
		frac_rvir = HOST_HALO_R200c/HOST_HALO_Rvir
		RADIUS = HOST_HALO_R200c
	if delta_crit == '500c' :
		frac_rvir = HOST_HALO_R500c/HOST_HALO_Rvir
		RADIUS = HOST_HALO_R500c
	if delta_crit == 'vir' :
		frac_rvir = HOST_HALO_Rvir/HOST_HALO_Rvir
		RADIUS = hdu_clu_bin['HOST_HALO_Rvir']
	if delta_crit == '2rvir' :
		frac_rvir = 2.*n.ones_like(HOST_HALO_R200c)
		RADIUS = 2.*hdu_clu_bin['HOST_HALO_Rvir']

	id_U = n.unique(hdu_clu['HOST_XRAY_image_path'], return_index = True )[1]
	id_U_130 = n.unique(hdu_clu['HOST_XRAY_image_path'][(HOST_HALO_M200c>=10**13.0) ], return_index = True )[1]
	id_U_135 = n.unique(hdu_clu['HOST_XRAY_image_path'][(HOST_HALO_M200c>=10**13.4) & (HOST_HALO_M200c<10**13.6)], return_index = True )[1]
	id_U_140 = n.unique(hdu_clu['HOST_XRAY_image_path'][(HOST_HALO_M200c>=10**14.0) & (HOST_HALO_M200c<10**14.2)], return_index = True )[1]
	id_U_145 = n.unique(hdu_clu['HOST_XRAY_image_path'][(HOST_HALO_M200c>=10**14.4) & (HOST_HALO_M200c<10**14.6)], return_index = True )[1]
	area_Mpc2 = n.sum(n.pi * (RADIUS[id_U] / 1000)**2)
	area_Mpc2_1e130 = n.sum(n.pi * (RADIUS[id_U_130] / 1000)**2)
	area_Mpc2_1e135 = n.sum(n.pi * (RADIUS[id_U_135] / 1000)**2)
	area_Mpc2_1e14 = n.sum(n.pi * (RADIUS[id_U_140] / 1000)**2)
	area_Mpc2_1e145 = n.sum(n.pi * (RADIUS[id_U_145] / 1000)**2)


	# SPIDERS DATA
	sel_spiders_host = (spiders_host['NMEMBERS'] >= 6) & (spiders_host['SCREEN_CLUZSPEC'] >= zmin) & (spiders_host['SCREEN_CLUZSPEC'] < zmax) & (spiders_host['M200c']>1e14)
	z_cluster_spiders = spiders_host['SCREEN_CLUZSPEC'][sel_spiders_host]
	z_cluster_spiders_mean = n.mean(spiders_host['SCREEN_CLUZSPEC'][sel_spiders_host])
	conversion = (1 /cosmo.arcsec_per_kpc_proper(z_cluster_spiders).to(u.deg /	u.Mpc)).value
	radius_Mpc_spiders = spiders_host['R200C_DEG'][sel_spiders_host] * conversion
	area_Mpc2_spiders = n.sum(n.pi * (radius_Mpc_spiders)**2)
	sel_spiders = (spiders_gal['ISMEMBER'] == 1) & (spiders_gal['SCREEN_CLUZSPEC'] >= zmin) & (spiders_gal['SCREEN_CLUZSPEC'] < zmax)
	abs_mag_i_spiders = spiders_gal['CMODELMAG_I'][sel_spiders] - cosmo.distmod(spiders_gal['SCREEN_CLUZSPEC'][sel_spiders]).value
	weight_spiders = 1./frac_obs_itp(spiders_gal['CMODELMAG_I'][sel_spiders])

	# theory input
	#total_number = n.cumsum(Schechter_M_z_M200c(mrange,z_cluster_mean,n.mean(HOST_HALO_M200c)) *area_Mpc2 *delta_mag)

	delta_mag = 0.01
	mrange = n.arange(-30, -10, delta_mag)
	mass_pc = scoreatpercentile(HOST_HALO_M200c, [25, 50, 75])
	LF_input_0 = Schechter_M_z_M200c(mrange,(zmin+zmax)/2., 10**13.5) #* n.log(10)#mass_pc[0])
	LF_input_1 = Schechter_M_z_M200c(mrange,(zmin+zmax)/2., 10**14.1) #* n.log(10)#mass_pc[1])
	LF_input_2 = Schechter_M_z_M200c(mrange,(zmin+zmax)/2., 10**14.5) #* n.log(10)#mass_pc[2])
	
	# figure of LF
	fig_out = os.path.join(fig_dir, 'spiders_cluster_galaxy_LF_z_'+aexp_str+'.png')
	mdex = 0.1
	mbins = n.arange(-30, -10, mdex)
	p.figure(1, (6., 5.5))
	p.axes([0.15, 0.15, 0.8, 0.8])
	if (zmin+zmax)/2.>=0.05 and (zmin+zmax)/2.<0.25:
		p.fill_between( mag_R18_015, y1 = phi_low_R18_015, y2 = phi_up_R18_015, alpha=0.2, label='R18 z=0.14' )
	if (zmin+zmax)/2.>=0.25 and (zmin+zmax)/2.<0.35:
		p.fill_between( mag_R18_030, y1 = phi_low_R18_030, y2 = phi_up_R18_030, alpha=0.2, label='R18 z=0.30' )
	p.plot(mrange, LF_input_0, 'b--', lw=1.5)#, label='R18 13.5') # +str(n.round(n.log10(mass_pc[0]),1)))
	p.plot(mrange, LF_input_1, 'k--', lw=1.5)#, label='R18 14.1') # +str(n.round(n.log10(mass_pc[1]),1)))
	p.plot(mrange, LF_input_2, 'r--', lw=1.5)#, label='R18 14.5') # +str(n.round(n.log10(mass_pc[2]),1)))
	# p.tight_layout()
	#p.hist(mag_abs_r[sel]    , lw=3, weights=n.ones(N_in_cluster)/mdex/area_Mpc2, bins=mbins, histtype='step', rasterized=True, label='magneticum, r')
	#s130 = (HOST_HALO_M200c>10**13.0) 
	#out = p.hist(mock_mag_abs[(x <= RADIUS)], lw=2, weights=n.ones_like(mock_mag_abs[(x <= RADIUS)]) / mdex / area_Mpc2 / n.log(10), bins = mbins, histtype = 'step', rasterized = True, label = r'All', color='green')
	#out = p.hist(mock_mag_abs[is_quiescent &(x <= RADIUS)], lw=2, weights=n.ones_like(mock_mag_abs[is_quiescent & (x <= RADIUS)]) / mdex / area_Mpc2 / n.log(10), bins = mbins, histtype = 'step', rasterized = True, label = r'All quiescent', color='green', ls='dotted')
	#out = p.hist(mock_mag_abs[x <= 0.00001], lw=2, weights=n.ones_like(mock_mag_abs[(x <= 0.00001)]) / mdex / area_Mpc2 / n.log(10), bins = mbins, histtype = 'step', rasterized = True, label = r'BCG', color='green', ls='dashed')
	#  HOST_HALO_M200c > 14.5
	#####
	s145 = (HOST_HALO_M200c>10**14.4) & (HOST_HALO_M200c<10**14.6)
	out = p.hist(mock_mag_abs[(x <= RADIUS)&(s145)], lw=2, 
			  weights=n.ones_like(mock_mag_abs[(x <= RADIUS)&(s145)]) / mdex / area_Mpc2_1e145 / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'$10^{14.4}<M_{200c}[M_\odot]<10^{14.6}$', color='red')
	out = p.hist(mock_mag_abs[is_quiescent & (x <= RADIUS)&(s145)], lw=2, 
			  weights=n.ones_like(mock_mag_abs[is_quiescent & (x <= RADIUS)&(s145)]) / mdex / area_Mpc2_1e145 / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'quiescent', color='red', ls='dotted')
	out = p.hist(mock_mag_abs[(x <= 0.0001)&(s145)], lw=2, 
			  weights=n.ones_like(mock_mag_abs[(x <= 0.0001)&(s145)]) / mdex / area_Mpc2_1e145 / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'BCG', color='red', ls='dashed')
	#####
	#  HOST_HALO_M200c > 14
	#####
	s14 = (HOST_HALO_M200c>10**14.0) & (HOST_HALO_M200c<10**14.2)
	out = p.hist(mock_mag_abs[(x <= RADIUS)&(s14)], lw=2, 
			  weights = n.ones_like(mock_mag_abs[(x <= RADIUS)&(s14)]) / mdex / area_Mpc2_1e14 / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'$10^{14}<M_{200c}[M_\odot]<10^{14.2}$', color='black')
	out = p.hist(mock_mag_abs[is_quiescent & (x <= RADIUS)&(s14)], lw=2, 
			  weights = n.ones_like(mock_mag_abs[is_quiescent & (x <= RADIUS)&(s14)]) / mdex / area_Mpc2_1e14 / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'quiescent', color='black', ls='dotted')
	out = p.hist(mock_mag_abs[(x <= 0.0001)&(s14)], lw=2, 
			  weights=n.ones_like(mock_mag_abs[(x <= 0.0001)&(s14)]) / mdex / area_Mpc2_1e14 / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'BCG', color='black', ls='dashed')	#####
	##
	# SPIDERS
	##
	dout = p.hist(abs_mag_i_spiders, lw=2, 
			  weights = weight_spiders / mdex / area_Mpc2_spiders / n.log(10), 
			  bins = mbins, histtype = 'step', rasterized = True, label = r'SPIDERS', color='green')

	p.title(r'$\bar{z}$=' + str(n.round(z_snap, 3)))
	p.ylabel(r'$\Phi=N_{gal} mag^{-1} Mpc^{-2}$')
	p.xlabel('Absolute magnitude r')
	p.yscale('log')
	p.xlim((n.min(abs_mag_i_spiders) - mdex,
			n.max(abs_mag_i_spiders) + mdex))
	p.ylim((1e-3, 20))
	p.legend(loc=2, frameon=False)
	p.grid()
	p.savefig(fig_out)
	p.clf()

	#fig_out = os.path.join(fig_dir, 'cluster_galaxy_z-mag-comparison.png')


	#p.figure(1, (6., 5.5))
	#p.tight_layout()
	#p.plot(zr_CLU, mock_mag_abs, 'k,', rasterized=True, label='HAM Ricci 18')
	#p.plot(z_spiders, abs_mag_i_spiders, 'b,', rasterized=True, label='SPIDERS')
	#p.legend(loc=0)
	## p.title(baseName)
	## p.yscale('log')
	#p.xlabel('redshift')
	#p.ylabel('mag')
	#p.grid()
	#p.savefig(fig_out)
	#p.clf()

