"""
What it does
------------

Computes magnitudes for galaxies in clusters following the luminosity function of Ricci et al. 2018

Adds color to red sequence galaxies
Model taken from redmapper trained on SDSS-IV/SPIDERS by J. Ider Chitham

Make a realization of galaxies following SDSS depth maps

outputs to a single file: $ENV_eRO_CLU_SAT_RS.fit

on ds43

References
----------

 * Ricci et al. 2018
 * Rykoff et al. 2013
 * Clerc et al. 2016

Command to run
--------------

python3 004_4_red_sequence.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

red sequence models

luminosity function

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
import numpy as n
from astropy.table import Table, Column
import kcor_lib as kk

print('Create file with galaxies around clusters')
print('=> Abundance matching for magnitudes')
print('=> Red sequence colors')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

env = sys.argv[1]
baseName = sys.argv[2]
delta_crit = sys.argv[3]
#env="MD10" 
#baseName="all_0.89510"
#delta_crit = '200c'
##delta_crit = '500c'
##delta_crit = 'vir'
##delta_crit = '2rvir'

z_snap = 1./float(baseName.split('_')[1])-1.
aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
print(env, baseName)
test_dir = os.path.join(os.environ[env])
path_2_CLU_SAT_catalog = os.path.join(test_dir, 'fits', baseName + '_galaxiesAroundClusters.fit')


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

fig_dir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'figures',
    env,
    'clusters',
)
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)

# ==============================
# READ red sequence DATA
# ==============================
red_seq_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'redSequence')
#path_2_f1 = os.path.join(red_seq_dir, 'sdss_cal_iter1_pars.fit')
#path_2_f2 = os.path.join(red_seq_dir, 'sdss_cal_zspec_redgals_model.fit')
path_2_f1 = os.path.join(red_seq_dir, 'sdss_cal_iter1_pars_2019.08.27.fit')
path_2_f2 = os.path.join(red_seq_dir, 'sdss_cal_zspec_redgals_model.2019.08.27.fit')

f1 = fits.open(path_2_f1)
f2 = fits.open(path_2_f2)

redshift_RS = f2[1].data['nodes'][0]
ug_RS = f2[1].data['meancol'].T[0].T[0]
gr_RS = f2[1].data['meancol'].T[1].T[0][:7]
ri_RS = f2[1].data['meancol'].T[2].T[0][:10]
iz_RS = f2[1].data['meancol'].T[3].T[0][:12]

ug_sigma_RS = f2[1].data['meancol_scatter'].T[0].T[0]
gr_sigma_RS = f2[1].data['meancol_scatter'].T[1].T[0][:7]
ri_sigma_RS = f2[1].data['meancol_scatter'].T[2].T[0][:10]
iz_sigma_RS = f2[1].data['meancol_scatter'].T[3].T[0][:12]

ug_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
                     n.hstack((ug_RS[0], ug_RS, ug_RS[-1])))
gr_RS_itp = interp1d(n.hstack((0., redshift_RS[:7], 2.0)),
                     n.hstack((gr_RS[0], gr_RS, gr_RS[-1])))
ri_RS_itp = interp1d(n.hstack((0., redshift_RS[:10], 2.0)),
                     n.hstack((ri_RS[0], ri_RS, ri_RS[-1])))
iz_RS_itp = interp1d(n.hstack((0., redshift_RS[:12], 2.0)),
                     n.hstack((iz_RS[0], iz_RS, iz_RS[-1])))

ug_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
    (ug_sigma_RS[0], ug_sigma_RS, ug_sigma_RS[-1])))
gr_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS[:7], 2.0)), n.hstack(
    (gr_sigma_RS[0], gr_sigma_RS, gr_sigma_RS[-1])))
ri_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS[:10], 2.0)), n.hstack(
    (ri_sigma_RS[0], ri_sigma_RS, ri_sigma_RS[-1])))
iz_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS[:12], 2.0)), n.hstack(
    (iz_sigma_RS[0], iz_sigma_RS, iz_sigma_RS[-1])))

"""
import matplotlib.pyplot as p

fig_out = os.path.join(red_seq_dir, 'red_sequence_2019_08_27.png')

p.figure(1, (6., 6.))

p.errorbar(redshift_RS[:7], gr_RS, xerr=0.05, yerr=gr_sigma_RS, label='g-r')
p.errorbar(redshift_RS[:10], ri_RS, xerr=0.05, yerr=ri_sigma_RS, label='r-i')
p.errorbar(redshift_RS[:12], iz_RS, xerr=0.05, yerr=iz_sigma_RS, label='i-z')

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
# ==============================
# Luminosity function model
# ==============================
# abundance matching for magnitudes based on Ricci et al. 2018. 
# Tables 3,5 
# Equation 8,9,10

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


# ==============================
# READ mock catalogs
# ==============================
# clusters
#hdu_host = fits.open(path_2_CLU_catalog)

# galaxies
hdu_clu = Table.read(path_2_CLU_SAT_catalog)  # , mode='update')
x = ((hdu_clu['x']-hdu_clu['HOST_HALO_x'])**2. + (hdu_clu['y']-hdu_clu['HOST_HALO_y'])**2. + (hdu_clu['z']-hdu_clu['HOST_HALO_z'])**2.)**0.5
is_quiescent = hdu_clu['is_quiescent']
zr_CLU = hdu_clu['redshift_R']
vmax = hdu_clu['vmax']
z_cluster = hdu_clu['HOST_redshift_R']
richness = n.zeros_like(zr_CLU)

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)

HOST_HALO_Mvir = hdu_clu['HOST_HALO_Mvir'] / h
HOST_HALO_Rvir = hdu_clu['HOST_HALO_Rvir']
HOST_HALO_M500c = hdu_clu['HOST_HALO_M500c'] / h
HOST_HALO_R500c = (DeltaVir_bn98(z_snap)/500. * HOST_HALO_M500c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir
HOST_HALO_M200c = hdu_clu['HOST_HALO_M200c'] / h
HOST_HALO_R200c = (DeltaVir_bn98(z_snap)/200. * HOST_HALO_M200c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir

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

N_gal = len(zr_CLU)

z_bins_2 = n.arange(0., n.max(zr_CLU) + 0.1, 0.001)
DM2_itp = interp1d(z_bins_2, cosmo.distmod(z_bins_2).value)

print('=============================================')
print('=============================================')
print('Abundance matching')
# ==============================
## ABUNDANCE MATCHING ##
# ==============================
# assigns magnitudes to the galaxies in the clusters
# loops as a function of halo mass 
# bins of d log 10 = 0.25
# 

log10M200c = n.log10(HOST_HALO_M200c)

dlog10M200c = 0.2
m200c_bins = n.arange( 12.6, 16.2, dlog10M200c )
for m200c_min, m200c_max in zip(m200c_bins, m200c_bins + dlog10M200c):
	print("======================================================")
	print('Mmin, Mmax', m200c_min, m200c_max)
	m_sel = (log10M200c>=m200c_min) & (log10M200c<m200c_max)
	N_gal = len(log10M200c[m_sel])
	print(m200c_min, m200c_max, N_gal)
	if N_gal>0:
		print('N galaxies', len(log10M200c[m_sel]))
		z_cluster_mean = n.mean(z_cluster[m_sel])
		print('z_cluster_mean', z_cluster_mean)
		print('N clusters=', len(n.unique(hdu_clu['HOST_HALO_id'][m_sel])))
		hdu_clu_bin = hdu_clu[m_sel]
		x = ((hdu_clu_bin['x']-hdu_clu_bin['HOST_HALO_x'])**2. + (hdu_clu_bin['y']-hdu_clu_bin['HOST_HALO_y'])**2. + (hdu_clu_bin['z']-hdu_clu_bin['HOST_HALO_z'])**2.)**0.5
		vmax = hdu_clu_bin['vmax']
		HOST_HALO_M200c = hdu_clu_bin['HOST_HALO_M200c'] / h
		print( len(n.unique(hdu_clu_bin['HOST_HALO_id'])) )
		print( len(n.unique(hdu_clu_bin['HOST_XRAY_image_path'])) )
		id_U = n.unique(hdu_clu_bin['HOST_XRAY_image_path'], return_index = True )[1]
		area_Mpc2 = n.sum(n.pi * (RADIUS[id_U] / 1000)**2)
		print('area_Mpc2', area_Mpc2)
		print('N galaxies around _cluster', N_gal)
		print('N galaxies witin X rvir cluster', len(HOST_HALO_M200c[x <= RADIUS[m_sel]]))
		# scatter in the abundance matching relation: random Gaussian number,
		# scale=0.15, loc=0
		rds = norm.rvs(loc=0, scale=0.15, size=N_gal)
		vmax_in_cluster = 10**(n.log10(vmax) + rds)
		delta_mag = 0.01
		mrange = n.arange(-30, -10, delta_mag)
		diff_number = Schechter_M_z_M200c(mrange, z_cluster_mean, n.mean(HOST_HALO_M200c)) * area_Mpc2 * delta_mag * n.log(10)
		total_number = n.cumsum(
			Schechter_M_z_M200c(
				mrange,
				z_cluster_mean,
				n.mean(HOST_HALO_M200c)) *
			area_Mpc2 * delta_mag * n.log(10))
		print('total_number', total_number, diff_number)
		x_itp = n.hstack((total_number[total_number > 0.1], total_number[total_number > 0.1][-1] * 10))
		y_itp = n.hstack((mrange[total_number > 0.1], mrange[total_number > 0.1][-1]))
		itp_clu = interp1d(x_itp, y_itp)
		mags_2_assign = itp_clu(n.arange(N_gal) + 1)
		print('mags_2_assign', mags_2_assign, mags_2_assign.shape)
		id_sort_vmax = n.argsort(vmax_in_cluster)[::-1]
		# id=0 the lowest vmax
		inter = n.zeros_like(vmax)
		#print('inter', inter)
		inter[id_sort_vmax] = mags_2_assign
		#print('inter', inter)
		mag_abs_r_HAM = inter
		#print('mag_abs_r_HAM', mag_abs_r_HAM)
		#print('dt=', time.time() - t0)
		print('HAM r magnitudes, time elapsed', mag_abs_r_HAM[:10], time.time() - t0)
		DM2 = DM2_itp(hdu_clu_bin['redshift_R'])
		hdu_clu['mag_abs_r'][m_sel] = mag_abs_r_HAM
		hdu_clu['mag_r'][m_sel] = mag_abs_r_HAM + DM2


x = ((hdu_clu['x']-hdu_clu['HOST_HALO_x'])**2. + (hdu_clu['y']-hdu_clu['HOST_HALO_y'])**2. + (hdu_clu['z']-hdu_clu['HOST_HALO_z'])**2.)**0.5
is_quiescent = hdu_clu['is_quiescent']
zr_CLU = hdu_clu['redshift_R']
#vmax = hdu_clu['vmax']
#z_cluster = hdu_clu['HOST_redshift_R']
richness = n.zeros_like(zr_CLU)

DM2 = DM2_itp(zr_CLU)

#omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
#DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)

#HOST_HALO_Mvir = hdu_clu['HOST_HALO_Mvir'] / h
#HOST_HALO_Rvir = hdu_clu['HOST_HALO_Rvir']
#HOST_HALO_M500c = hdu_clu['HOST_HALO_M500c'] / h
#HOST_HALO_R500c = (DeltaVir_bn98(z_snap)/500. * HOST_HALO_M500c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir
#HOST_HALO_M200c = hdu_clu['HOST_HALO_M200c'] / h
#HOST_HALO_R200c = (DeltaVir_bn98(z_snap)/200. * HOST_HALO_M200c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir

N_gal = len(zr_CLU)
# ==============================
## COMPUTES THE RICHNESS 
# ==============================
# select galaxies brigther than M star + 2, -20
m_star_p1 = -22 + 2.0
bright_red_galaxies = (hdu_clu['mag_abs_r'] <= m_star_p1) & (is_quiescent) & (x <= RADIUS)
faint_galaxies = (~bright_red_galaxies)
test3 = n.unique(
    hdu_clu['HOST_HALO_id'][bright_red_galaxies],
    return_counts=True,
    return_index=True,
    return_inverse=True)
richness_per_cluster = interp1d(test3[0], test3[3])
richness_i = richness_per_cluster(hdu_clu['HOST_HALO_id'][bright_red_galaxies])
richness[bright_red_galaxies] = richness_i
richness_col = Column(
    name='richness',
    dtype=n.int64,
    unit='',
    data = richness)
print('Richness estimate (M star + 2): ', richness[:10], time.time() - t0)

# ==============================
## COMPUTES THE COLOR ##
# ==============================
#rds = norm.rvs(loc=0, scale=0.05, size=len(zr_CLU[is_quiescent]))
rds = norm.rvs(loc=0, scale=1., size=len(zr_CLU[is_quiescent]))

scatter_gr = gr_sigma_RS_itp (zr_CLU[is_quiescent])
gr = gr_RS_itp(zr_CLU[is_quiescent]) + rds * scatter_gr

scatter_ri = ri_sigma_RS_itp (zr_CLU[is_quiescent])
ri = ri_RS_itp(zr_CLU[is_quiescent]) + rds * scatter_ri

scatter_iz = iz_sigma_RS_itp (zr_CLU[is_quiescent]) 
iz = iz_RS_itp(zr_CLU[is_quiescent]) + rds * scatter_iz

gr_all = n.zeros_like(zr_CLU)
ri_all = n.zeros_like(zr_CLU)
iz_all = n.zeros_like(zr_CLU)
gr_all[is_quiescent] = gr
ri_all[is_quiescent] = ri
iz_all[is_quiescent] = iz

gr_col = Column(name='galaxy_gr', dtype=n.float32, unit='mag', data = gr_all)
ri_col = Column(name='galaxy_ri', dtype=n.float32, unit='mag', data = ri_all)
iz_col = Column(name='galaxy_iz', dtype=n.float32, unit='mag', data = iz_all)
print('colors assigned',
      gr_all[:10],
      ri_all[:10],
      iz_all[:10],
      time.time() - t0)

# ==============================
## COMPUTES SDSS magnitudes ##
# ==============================
# Add uncertainties on the magnitudes and colors following SDSS depth maps
# reads the maps
sdss_depth_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'sdss_depth')
depth_sdss_g = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_g_model_10sigma.fits.gz')
depth_sdss_i = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_i_model_10sigma.fits.gz')
depth_sdss_r = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_r_model_10sigma.fits.gz')
depth_sdss_z = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_z_model_10sigma.fits.gz')
#
NSIDE = 2048
sdss_g_depth = hp.fitsfunc.read_map(depth_sdss_g)
sdss_r_depth = hp.fitsfunc.read_map(depth_sdss_r)
sdss_i_depth = hp.fitsfunc.read_map(depth_sdss_i)
sdss_z_depth = hp.fitsfunc.read_map(depth_sdss_z)
# retrieves the depths
pixels = healpy.ang2pix(
    NSIDE,
    hdu_clu['RA'],
    hdu_clu['DEC'],
    nest=False,
    lonlat=True)
depth_10_sigma_g = sdss_g_depth[pixels]
depth_10_sigma_r = sdss_r_depth[pixels]
depth_10_sigma_i = sdss_i_depth[pixels]
depth_10_sigma_z = sdss_z_depth[pixels]
# computes the magnitudes
# K correction

#kcorr_value_g = kk.calc_kcor('g', zr_CLU, 'g - r', gr_all)
kcorr_value_r = kk.calc_kcor('r', zr_CLU, 'g - r', gr_all)
#kcorr_value_i = kk.calc_kcor('i', zr_CLU, 'g - i', gr_all+ri_all)
#kcorr_value_z = kk.calc_kcor('z', zr_CLU, 'g - z', gr_all+ri_all+iz_all)
# 
sdss_r_temp = hdu_clu['mag_abs_r'] + DM2 + kcorr_value_r
sdss_g_temp = sdss_r_temp + gr_all 
sdss_i_temp = sdss_r_temp - ri_all 
sdss_z_temp = sdss_r_temp - ri_all - iz_all 
# conver to flux
def get_sigma_m(depth_10_sigma, sdss_temp):
	limiting_flux = 2.5 * 10**(-0.4 * (depth_10_sigma)) / n.log(10)/10.
	observed_flux = 2.5 * 10**(-0.4 * (sdss_temp)) / n.log(10)
	return 2.5 * limiting_flux / observed_flux / n.log(10)

# their uncertainty
sdss_g_err = get_sigma_m(depth_10_sigma_g, sdss_g_temp)
sdss_r_err = get_sigma_m(depth_10_sigma_r, sdss_r_temp)
sdss_i_err = get_sigma_m(depth_10_sigma_i, sdss_i_temp)
sdss_z_err = get_sigma_m(depth_10_sigma_z, sdss_z_temp)
# assign uncertainty=0 outside of the SDSS footprint
sdss_g_err[sdss_g_err == n.inf] = 0
sdss_r_err[sdss_r_err == n.inf] = 0
sdss_i_err[sdss_i_err == n.inf] = 0
sdss_z_err[sdss_z_err == n.inf] = 0
print('SDSS magnitude uncertainties computed',
      sdss_g_err[:10],
      sdss_r_err[:10],
      sdss_i_err[:10],
      sdss_z_err[:10],
      time.time() - t0)

#s1 = (sdss_i_temp>21.3)&(sdss_i_temp<21.4)&(sdss_i_err>0)
#print(n.median(sdss_i_err[s1]) )

# random variable in a norm distribution with sigma = 1

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_g_col = Column(
    name='sdss_g',
    dtype=n.float32,
    unit='mag',
    data = sdss_g_temp +
    rds *
    sdss_g_err)

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_r_col = Column(
    name='sdss_r',
    dtype=n.float32,
    unit='mag',
    data = sdss_r_temp +
    rds *
    sdss_r_err)

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_i_col = Column(
    name='sdss_i',
    dtype=n.float32,
    unit='mag',
    data = sdss_i_temp +
    rds *
    sdss_i_err)

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_z_col = Column(
    name='sdss_z',
    dtype=n.float32,
    unit='mag',
    data = sdss_z_temp +
    rds *
    sdss_z_err)


sdss_g_err_col = Column(
    name='sdss_g_err',
    dtype=n.float32,
    unit='mag',
    data = sdss_g_err)
sdss_r_err_col = Column(
    name='sdss_r_err',
    dtype=n.float32,
    unit='mag',
    data = sdss_r_err)
sdss_i_err_col = Column(
    name='sdss_i_err',
    dtype=n.float32,
    unit='mag',
    data = sdss_i_err)
sdss_z_err_col = Column(
    name='sdss_z_err',
    dtype=n.float32,
    unit='mag',
    data = sdss_z_err)


print('creates hdu, time elapsed', time.time() - t0)


hdu_clu.add_column(gr_col)
hdu_clu.add_column(ri_col)
hdu_clu.add_column(iz_col)
hdu_clu.add_column(richness_col)
#hdu_clu.add_column(rmag_abs_ham_col)
#hdu_clu.add_column(rmag_ham_col)

hdu_clu.add_column(sdss_g_col)
hdu_clu.add_column(sdss_r_col)
hdu_clu.add_column(sdss_i_col)
hdu_clu.add_column(sdss_z_col)

hdu_clu.add_column(sdss_g_err_col)
hdu_clu.add_column(sdss_r_err_col)
hdu_clu.add_column(sdss_i_err_col)
hdu_clu.add_column(sdss_z_err_col)

hdu_clu.write(path_2_CLU_SAT_catalog, overwrite=True)

