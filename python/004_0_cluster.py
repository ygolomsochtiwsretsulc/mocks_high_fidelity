"""
What it does 
------------

Computes the X-ray cluster model in each shell

Uses scaling relations to infer LX_soft_cin, LX_soft_cex, TX_cin, TX_cex, FX_soft_cin, FX_soft_cex,  

References
----------

* Bulbul et al. 2019 https://ui.adsabs.harvard.edu/abs/2019ApJ...871...50B

Command to run
--------------

python3 004_0_cluster.py environmentVAR fileBasename

arguments 
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, extinction, matplotlib

"""
print('Creates the h5 Cluster file per shell')
print('------------------------------------------------')
print('------------------------------------------------')
import time 
t0 = time.time()
import numpy as n
import os, sys
import h5py
from astropy_healpix import healpy 

from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.stats import norm
from scipy.special import erf

import astropy.io.fits as fits
import astropy.constants as cc
import astropy.units as u

# import all pathes 
import os, sys

env = sys.argv[1] # 'MD04'
baseName = sys.argv[2] # "all_0.62840"
print(env, baseName)
make_figure = True
make_figure = False

# initializes pathes to files
test_dir = os.path.join( os.environ[env], 'hlists', 'fits' )

path_2_light_cone = os.path.join(test_dir, baseName+'.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName+'_coordinates.h5')
path_2_galaxy_file = os.path.join(test_dir, baseName+'_galaxy.h5')
path_2_CLU_file = os.path.join(test_dir, baseName+'_CLU.h5')

# x ray extinction for clusters, from A. Finoguenov
path_2_NH_law = os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", 'nh_flux.tbl' )
NH_DATA = n.loadtxt(path_2_NH_law , unpack=True)
nh_law = interp1d(n.hstack((-1.e30, 1.0*10**10, NH_DATA[0],1.0*10**30)), n.hstack((1., 1., NH_DATA[1][0]/NH_DATA[1],  7.03612982e-09 )))

# eRosita flux limits for point sources
path_2_flux_limits = os.path.join( os.environ['GIT_AGN_MOCK'], "data", "erosita", "flux_limits.fits")

# simulation setup 
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
if env=="MD10" or env=="MD04":
	cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)
	h = 0.6777
	L_box = 1000.0 / h
	cosmo = cosmoMD
if env=="UNIT_fA1_DIR" or env=="UNIT_fA1i_DIR" or env=="UNIT_fA2_DIR":
	cosmoUNIT = FlatLambdaCDM(H0=67.74*u.km/u.s/u.Mpc, Om0=0.308900)
	h = 0.6774
	L_box = 1000.0 / h
	cosmo = cosmoUNIT

f1 = fits.open(path_2_light_cone)
Mvir = f1[1].data['Mvir']/h 
M500c = f1[1].data['M500c']/h 
log_vmax = n.log10(f1[1].data['vmax'])
N_obj = len(Mvir)
cluster = ( M500c > 5e13 )
scale_of_last_MM = f1[1].data['scale_of_last_MM'][cluster]
N_clu = len(Mvir[cluster])
f1.close()

ids_cluster = n.arange(N_obj)[cluster]

f2 = h5py.File(path_2_coordinate_file, 'r')
ra = f2['/coordinates/ra'].value            [cluster]
dec = f2['/coordinates/dec'].value          [cluster]
zz = f2['/coordinates/redshift_R'].value    [cluster]
zzs = f2['/coordinates/redshift_S'].value   [cluster]
dL_cm = f2['/coordinates/dL'].value         [cluster]
galactic_NH = f2['/coordinates/NH'].value   [cluster]
galactic_ebv = f2['/coordinates/ebv'].value [cluster]
g_lat = f2['/coordinates/g_lat'].value      [cluster]
g_lon = f2['/coordinates/g_lon'].value      [cluster]
ecl_lat = f2['/coordinates/ecl_lat'].value  [cluster]
ecl_lon = f2['/coordinates/ecl_lon'].value  [cluster]
f2.close()
print('coordinate file opened', time.time()-t0)

# cosmological volume
zmin = n.min(zz)
zmax = n.max(zz)
z_mean = 0.5*(zmin+zmax)
print(zmin, '<z<', zmax)
vol = (cosmo.comoving_volume(zmax).value-cosmo.comoving_volume(zmin).value)
DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
print('volume',vol,'Mpc3', time.time()-t0)

f3 = h5py.File(path_2_galaxy_file, 'r')
mass = f3['/galaxy/SMHMR_mass'].value # log of the stellar mass
#f3.close()


# Bulbul et al. 2018 scaling relations
z_pivot = 0.45
M_pivot = 6.35*10**(14)
E_pivot = cosmo.efunc(z_pivot)
#
scatter_1 = norm.rvs(loc=0, scale=1., size=N_clu)
#
# Eq 17 T X cin
SR_0 = lambda m500, z : 6.41 * (m500/M_pivot) ** 0.799 * (cosmo.efunc(z)/E_pivot)**(2./3.) * ((1+z)/(1+z_pivot))**(-0.36)
# Eq 18 T X cex
SR_1 = lambda m500, z : 6.09 * (m500/M_pivot) ** 0.799 * (cosmo.efunc(z)/E_pivot)**(2./3.) * ((1+z)/(1+z_pivot))**(-0.33)
# Eq 20 M ICM cex
SR_2 = lambda m500, z : 7.37 * 1e13 * (m500/M_pivot) ** 1.259 * ((1+z)/(1+z_pivot))**(0.18)
# Eq 25 L X cin
SR_3 = lambda m500, z : 4.12 * 1e44 * (m500/M_pivot) ** 1.89 * (cosmo.efunc(z)/E_pivot)**2.0 * ((1+z)/(1+z_pivot))**(-0.2)
# Eq 26 L X cex
SR_4 = lambda m500, z : 2.84 * 1e44 * (m500/M_pivot) ** 1.60 * (cosmo.efunc(z)/E_pivot)**2.0 * ((1+z)/(1+z_pivot))**(-0.1)
#
TX_BB_18_cin = SR_0(M500c[cluster], zz)
TX_BB_18_cex = SR_1(M500c[cluster], zz)
MICM_BB_18_cex = SR_2(M500c[cluster], zz)
LX_BB_18_cin = SR_3(M500c[cluster], zz)
LX_BB_18_cex = SR_4(M500c[cluster], zz)
print('LX_BB_18_cin', LX_BB_18_cin[:20], time.time()-t0)
print('LX_BB_18_cex', LX_BB_18_cex[:20], time.time()-t0)

# Flux limit for point sources
# use the flux limit for SNR3 point sources to have all possible clusters, in particular at high redshift
# the a second flux limit will be applied with image simulations
pix_ids = healpy.ang2pix(512, n.pi/2.-g_lat*n.pi/180., g_lon*n.pi/180., nest=True)
flux_lim_data = fits.open(path_2_flux_limits)
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3  
flux_limit = flux_lim_data[1].data['flux_limit_SNR3'][pix_ids]
print('flux limit file opened', time.time()-t0)

# Flux observed and attenuated 
attenuation = nh_law(galactic_NH)
CLU_FX_soft = LX_BB_18_cin / (4*n.pi*dL_cm**2.)
CLU_FX_soft_attenuated = attenuation * CLU_FX_soft
detected = (CLU_FX_soft_attenuated>10**(flux_limit))

# relaxation state
# the ones that have the most recent MM are most disturbed.
# we do 4 classes with equal number of each class
time_of_last_MM = cosmo.age(1/scale_of_last_MM-1).value
time_now = cosmo.age(zz).value
delta_t_MM = time_now - time_of_last_MM
#
delta_t_hist = n.histogram(delta_t_MM, bins=100)
delta_t_fraction = n.hstack(( 0., n.cumsum(delta_t_hist[0])/N_clu ))
delta_t_values_itp = interp1d( delta_t_hist[1], delta_t_fraction)
# if coolness is small, then it is disturbed
# if coolness is long, then it is relaxed
coolness = delta_t_values_itp(delta_t_MM)


# writes the results
print('writes results', time.time()-t0)
f = h5py.File(path_2_CLU_file, "a")
f.attrs['file_name'] = os.path.basename(path_2_CLU_file)
f.attrs['creator']   = 'JC'

# writes the results
halo_data = f.create_group( 'CLUSTERS' )

halo_data.create_dataset( 'ids_cluster', data = ids_cluster )

ds = halo_data.create_dataset( 'LX_soft_cin', data = n.log10(LX_BB_18_cin) + scatter_1 * 0.27 )
ds.attrs['units'] = 'log10(L_X/[0.5-2keV, cin, erg/s])'

halo_data.create_dataset('LX_soft_cex', data = n.log10(LX_BB_18_cex) + scatter_1 * 0.27 )
ds.attrs['units'] = 'log10(L_X/[0.5-2keV, cex, erg/s])'

halo_data.create_dataset('TX_cin', data = n.log10(TX_BB_18_cin) + scatter_1 * 0.179 )
ds.attrs['units'] = 'log10(T cin [keV])'

halo_data.create_dataset('TX_cex', data = n.log10(TX_BB_18_cex) + scatter_1 * 0.128 )
ds.attrs['units'] = 'log10(T cex [keV])'

halo_data.create_dataset('MICM_cex', data = n.log10(MICM_BB_18_cex) + scatter_1 * 0.098 )
ds.attrs['units'] = 'log10(M_ICM [Msun]))'

halo_data.create_dataset('FX_soft', data = CLU_FX_soft )
ds.attrs['units'] = 'F_X / [0.5-2keV, erg/cm2/s]'

halo_data.create_dataset('FX_soft_attenuated', data = CLU_FX_soft_attenuated )
ds.attrs['units'] = 'F_X / [0.5-2keV, erg/cm2/s]'

halo_data.create_dataset('detected', data = detected )

halo_data.create_dataset('scatter_1', data = scatter_1 )

halo_data.create_dataset('coolness', data = coolness )
ds.attrs['units'] = '0:very disturbed, 1:relaxed'

f.close()

### NOW THE FIGURES ###
if make_figure:
	import matplotlib
	matplotlib.use('Agg')
	matplotlib.rcParams.update({'font.size': 14})
	import matplotlib.pyplot as p


	fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'clusters', )
	if os.path.isdir(fig_dir)==False:
		os.system('mkdir -p '+fig_dir)

	fig_out = os.path.join(fig_dir, 'Cluster_LX_soft_cin_vs_mass_'+baseName+'.png')
	X = M500c[cluster]
	Y = n.log10(LX_BB_18_cin) + scatter_1 * 0.27
	p.figure(1, (6.,5.5))
	p.tight_layout()
	p.plot(X, Y, 'k+', rasterized=True)
	p.title(baseName)
	p.xscale('log')
	p.xlabel('M500c')
	p.ylabel('LX_soft cin')
	p.grid()
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(fig_dir, 'Cluster_LX_soft_cex_vs_mass_'+baseName+'.png')
	X = M500c[cluster]
	Y = n.log10(LX_BB_18_cex) + scatter_1 * 0.27
	p.figure(1, (6.,5.5))
	p.tight_layout()
	p.plot(X, Y, 'k+', rasterized=True)
	p.title(baseName)
	p.xscale('log')
	p.xlabel('M500c')
	p.ylabel('LX_soft cex')
	p.grid()
	p.savefig(fig_out)
	p.clf()

sys.exit()

# other scaling relations

# DIETRICH 2018
N_Mgas = 31.92 # Dietrich 18
N_kT 	= 2.18
N_L 	= 103.7
N_Lce 	= 102.66

slope_E_Mgas 	= 0.05 # Dietrich 18
slope_E_kT 	= 0.61
slope_E_L 	= 1.20
slope_E_Lce 	= 1.82

slope_M500_Mgas= 1.398 # Dietrich 18
slope_M500_kT 	= 0.66
slope_M500_L 	= 1.43 # 1.26*(1.+0.33*0.43)
slope_M500_Lce = 1.36 # 1.06*(1.+0.33*0.88)

scatter_Mgas = 0.106 # Dietrich 18
scatter_kT = 0.18
scatter_L = 0.24
scatter_Lce = 0.17

# MANTZ 2016

#N_Mgas = 31.98 
#N_kT 	= 2.18
#N_L 	= 103.7
#N_Lce 	= 102.66

#slope_E_Mgas 	= -0.11 
#slope_E_kT 	= 0.61
#slope_E_L 	= 1.20
#slope_E_Lce 	= 1.82

#slope_M500_Mgas= 1.04
#slope_M500_kT 	= 0.66
#slope_M500_L 	= 1.26
#slope_M500_Lce = 1.06

#scatter_Mgas = 0.086
#scatter_kT = 0.18
#scatter_L = 0.24
#scatter_Lce = 0.17

E035 = cosmoMD.efunc(0.35)

# converts logM500 to clusters observables
m500_to_qty = lambda logM500, z, slope_efunc, slope_m500, normalization : n.e**normalization * (cosmoMD.efunc(z)/E035)**(slope_efunc) * (10**(logM500-n.log10(6)-14))**(slope_m500)

logM500_to_logMgas 	= lambda logM500, z : m500_to_qty( logM500, z, slope_E_Mgas, slope_M500_Mgas, N_Mgas)
logM500_to_kT 		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_kT, slope_M500_kT, N_kT)
logM500_to_L 		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_L, slope_M500_L, N_L)
logM500_to_Lce		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_Lce, slope_M500_Lce, N_Lce)

z = f1.attrs['redshift']
m500c_i = f1['/halo_properties/M500c'].value
ok = (m500c_i>1e13)
log_m500c = n.log10(m500c_i[ok])

nCluster = len(log_m500c)
#rds = (n.random.rand(len(log_m500c))-0.5)*2.                   

Mean_Mgas = n.log10(logM500_to_logMgas	(log_m500c, z))
V_scatter_Mgas = norm.rvs(loc=0,scale=scatter_Mgas,size=nCluster)
VAL_Mgas = n.zeros_like(m500c_i)
VAL_Mgas[ok] = Mean_Mgas + V_scatter_Mgas

Mean_kT = logM500_to_kT(log_m500c, z)
V_scatter_kT = norm.rvs(loc=0,scale=scatter_kT,size=nCluster)
VAL_kT= n.zeros_like(m500c_i)
VAL_kT[ok] = Mean_kT + V_scatter_kT

Mean_L = n.log10(logM500_to_L(log_m500c, z))
V_scatter_L = norm.rvs(loc=0,scale=scatter_L,size=nCluster)
VAL_L= n.zeros_like(m500c_i)
VAL_L[ok] = Mean_L + V_scatter_L

Mean_Lce = n.log10(logM500_to_Lce(log_m500c, z))
V_scatter_Lce = norm.rvs(loc=0,scale=scatter_Lce,size=nCluster)
VAL_Lce= n.zeros_like(m500c_i)
VAL_Lce[ok] = Mean_Lce + V_scatter_Lce
