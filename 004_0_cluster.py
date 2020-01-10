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


244 objects : 

182 XXL, 
49 HIFLUGS
12 XCOP


# correlated random variables 
# https://scipy-cookbook.readthedocs.io/items/CorrelatedRandomSamples.html
# import numpy as np
from scipy.linalg import eigh, cholesky
from scipy.stats import norm
# Choice of cholesky or eigenvector method.
method = 'cholesky'
#method = 'eigenvectors'
# The desired covariance matrix.
r = np.array([
        [  0.36, 0.49],
        [  0.49,  0.19],
    ])

# Generate samples from three independent normally distributed random
# variables (with mean 0 and std. dev. 1).
x = norm.rvs(size=(3, N_clu))

# We need a matrix `c` for which `c*c^T = r`.  We can use, for example,
# the Cholesky decomposition, or the we can construct `c` from the
# eigenvectors and eigenvalues.

if method == 'cholesky':
    # Compute the Cholesky decomposition.
    c = cholesky(r, lower=True)

# Convert the data to correlated random variables. 
y = np.dot(c, x)



"""
from astropy.table import Table, Column
from astropy_healpix import healpy
import sys
import os
import time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import h5py
import numpy as n
from lib_cluster import write_img, create_matrix
print('Creates the h5 Cluster file per shell')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


# import all pathes

env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "all_0.62840"
z_snap = 1./float(baseName.split('_')[1])-1.
print(env, baseName,z_snap)
make_figure = True
make_figure = False
is_writing_images = sys.argv[3] # 'yes'

# initializes pathes to files
test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_CLU_file = os.path.join(test_dir, baseName + '_CLU.fits')

# x ray extinction for clusters, from A. Finoguenov
path_2_NH_law = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "xray_k_correction",
    'nh_flux.tbl')
NH_DATA = n.loadtxt(path_2_NH_law, unpack=True)
nh_law = interp1d(
    n.hstack(
        (-1.e30,
         1.0 * 10**10,
         NH_DATA[0],
         1.0 * 10**30)),
    n.hstack(
            (1.,
             1.,
             NH_DATA[1][0] / NH_DATA[1],
             7.03612982e-09)))

# eRosita flux limits for point sources
path_2_flux_limits = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "erosita",
    "flux_limits.fits")

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

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)


f1 = fits.open(path_2_light_cone)
N_obj = len(f1[1].data['M500c'])
cluster = ( n.log10(f1[1].data['M500c'] / h) > 13 ) & ( f1[1].data['pid'] == -1 )

Mvir = f1[1].data['Mvir'][cluster] / h
Rvir = f1[1].data['Rvir'][cluster]
M500c = f1[1].data['M500c'][cluster] / h
R500c = (DeltaVir_bn98(z_snap)/500. * M500c / Mvir)**(1./3.)*Rvir
logM500c = n.log10(M500c)
log_vmax = n.log10(f1[1].data['vmax'][cluster])
scale_of_last_MM = f1[1].data['scale_of_last_MM'][cluster]
N_clu = len(Mvir)
print(N_clu, 'clusters')
Xoff = f1[1].data['Xoff'][cluster] / f1[1].data['Rvir'][cluster]
b_to_a_500c = f1[1].data['b_to_a_500c'][cluster]
#f1.close()
halo_lineID = n.arange(N_obj)[cluster]
ids_cluster = (n.arange(N_obj)[cluster]*1e12 + f1[1].data['id'][cluster] ).astype('int')
ids_cluster_str = n.array([str(el).zfill(22) for el in ids_cluster])

f2 = fits.open(path_2_coordinate_file)
ra = f2[1].data['ra'][cluster]
dec = f2[1].data['dec'][cluster]
zz = f2[1].data['redshift_R'][cluster]
zzs = f2[1].data['redshift_S'][cluster]
dL_cm = f2[1].data['dL'][cluster]
galactic_NH = f2[1].data['nH'][cluster]
galactic_ebv = f2[1].data['ebv'][cluster]
g_lat = f2[1].data['g_lat'][cluster]
g_lon = f2[1].data['g_lon'][cluster]
ecl_lat = f2[1].data['ecl_lat'][cluster]
ecl_lon = f2[1].data['ecl_lon'][cluster]
#f2.close()
print('coordinate file opened', time.time() - t0)

# cosmological volume
zmin = n.min(zz)
zmax = n.max(zz)
z_mean = 0.5 * (zmin + zmax)
print(zmin, '<z<', zmax)
vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value)
DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
print('volume', vol, 'Mpc3', time.time() - t0)

f3 = fits.open(path_2_galaxy_file)
mass = f3[1].data['SMHMR_mass'][cluster]  # log of the stellar mass
# f3.close()

path_2_cbp = os.path.join(os.environ['GIT_CBP'])

Mpc=3.0856776e+24
msun=1.98892e33
if N_clu < 100:
	nsim = N_clu*100
else:
	nsim = N_clu*10
covor=n.loadtxt(os.path.join(path_2_cbp, 'covmat_xxl_hiflugcs_xcop.txt'))
xgrid_ext=n.loadtxt(os.path.join(path_2_cbp, 'radial_binning.txt'))
mean_log=n.loadtxt(os.path.join(path_2_cbp, 'mean_pars.txt'))
coolfunc=n.loadtxt(os.path.join(path_2_cbp, 'coolfunc.dat'))

def calc_lx(prof,kt,m5,z):
	"""
	Compute the X-ray luminosity in the profile
	to be extended to 3x r500c
	"""
	ez2=cosmo.efunc(z)**2
	rhoc = cosmo.critical_density(z).value
	r500 = n.power(m5*msun/4.*3./n.pi/500./rhoc,1./3.)
	resfact=n.sqrt(kt/10.0)*n.power(ez2,3./2.)
	prof_em=prof*resfact # emission integral
	tlambda=n.interp(kt,coolfunc[:,0],coolfunc[:,1]) # cooling function
	dx=n.empty(len(xgrid_ext))
	dx[0]=xgrid_ext[0]
	dx[1:len(xgrid_ext)]=(n.roll(xgrid_ext,-1)-xgrid_ext)[:len(xgrid_ext)-1]
	#print(prof_em*xgrid_ext*r500**2*2.*n.pi*tlambda*Mpc*dx)
	lxcum=n.cumsum(prof_em*xgrid_ext*r500**2*2.*n.pi*tlambda*Mpc*dx) # riemann integral
	lx_500=n.interp(1.,xgrid_ext,lxcum) # evaluated at R500
	return lx_500

profs=n.exp(n.random.multivariate_normal(mean_log,covor,size=nsim))

allz=profs[:,len(mean_log)-3]
allkt=profs[:,len(mean_log)-1]
allm5=profs[:,len(mean_log)-2]

profiles=profs[:,:len(xgrid_ext)]

alllx=n.empty(nsim)
for i in range(nsim):
	tprof=profiles[i]
	alllx[i]=calc_lx(tprof,allkt[i],allm5[i],allz[i])


print('profiles created', time.time()-t0)
from sklearn.neighbors import BallTree

Tree_profiles = BallTree(n.transpose([allz, n.log10(allm5)]))
#DATA = n.array([[z_i, m_i] for z_i, m_i in zip(zzs, logM500c) ])
DATA = n.transpose([zzs, logM500c])
ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)

ids = n.hstack((ids_out))

# LX_out 0.5-2 rest frame erg/s
# to convert to 0.5-2 observed frame 
LX_out = alllx[ids]
KT_OUT = allkt[ids]
print('profiles matched', time.time()-t0)


dir_2_SMPT_image = os.path.join(os.environ[env], "cat_CLU_SIMPUT", 'cluster_images')
if os.path.isdir(dir_2_SMPT_image) == False:
    os.system('mkdir -p ' + dir_2_SMPT_image)

print('creates SIMPUT images')
path_2_images = []
# writes images using these profiles

for profile_i, r500c_i, conversion_arcmin, file_name, b_a in zip(profiles[ids], R500c,  cosmo.kpc_proper_per_arcmin(zz).value, ids_cluster_str, b_to_a_500c):
	image_file = os.path.join(dir_2_SMPT_image, file_name+'.fits')
	print(image_file)
	path_2_images.append(file_name)
	x_coord = n.hstack(( 0., xgrid_ext*r500c_i/conversion_arcmin, 10*xgrid_ext[-1]*r500c_i/conversion_arcmin, 1000*xgrid_ext[-1]*r500c_i/conversion_arcmin )) # arcminutes
	y_coord = n.hstack(( profile_i[0], profile_i, 0., 0. ))
	profile = interp1d(x_coord, y_coord)
	matrix = create_matrix(profile, n_pixel = 120, b_a = b_a)
	if is_writing_images == 'yes':
		write_img(matrix, image_file, n_pixel = 120)


# attenuation grid should cover down to 0.05 keV 
itp_logNH, itp_z, itp_kt, itp_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters.txt"), unpack=True )

nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
z_vals = n.hstack(( n.arange(0.,0.7,0.05), [0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6] ))
#kT_vals = n.hstack(( n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))
kT_vals = n.hstack(([0.1, 0.2], n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))

XX_nh, YY_z, ZZ_kt = n.meshgrid(nh_vals, z_vals, kT_vals)

shape_i = XX_nh.shape

matrix = itp_frac_obs.reshape(shape_i)

from scipy.interpolate import RegularGridInterpolator
attenuation_3d = RegularGridInterpolator((z_vals, n.log10(nh_vals*1e22), kT_vals), matrix)

sel = (galactic_NH>1e19)
bad_NH = (sel==False)
galactic_NH[bad_NH] = 1.000001e19

k_correction_3d = attenuation_3d( n.transpose([zzs, n.ones_like(zzs)*20.1, KT_OUT]))

itp_z, itp_kt, itp_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters_no_nH.txt"), unpack=True )

YY_z, ZZ_kt = n.meshgrid(z_vals, kT_vals)

shape_i = YY_z.shape

matrix_2d = itp_frac_obs.reshape(shape_i)

from scipy.interpolate import RegularGridInterpolator
attenuation_2d = RegularGridInterpolator((kT_vals, z_vals), matrix_2d)

k_correction_2d = attenuation_2d( n.transpose([KT_OUT, zzs]))


attenuate_X_logNH, attenuate_Y_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "nh_attenuation_clusters.txt"), unpack=True )

#att = attenuation( n.transpose([zzs, n.log10(data['galactic_NH']), KT_OUT]))

LX_obsF = LX_out / k_correction_2d 
CLU_FX_soft_1 = LX_obsF / (4 * n.pi * dL_cm**2.)

itp_attenuation = interp1d(attenuate_X_logNH, attenuate_Y_frac_obs)
att = itp_attenuation(n.log10(galactic_NH))
CLU_FX_soft = CLU_FX_soft_1 * att

# Flux limit for point sources
# use the flux limit for SNR3 point sources to have all possible clusters, in particular at high redshift
# the a second flux limit will be applied with image simulations
pix_ids = healpy.ang2pix(
    512,
    n.pi /
    2. -
    g_lat *
    n.pi /
    180.,
    g_lon *
    n.pi /
    180.,
    nest=True)
flux_lim_data = fits.open(path_2_flux_limits)
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3
flux_limit = flux_lim_data[1].data['flux_limit_SNR3'][pix_ids]
print('flux limit file opened', time.time() - t0)
detected = (CLU_FX_soft > 10**(flux_limit))

# relaxation state
# based on MD40, quartiles of Xoff
xoff_bins_Q = n.array([0.,  4.28608992e-02, 6.59490117e-02, 1.00837135e-01, 5.7516e-01])
coolness = n.ones_like(Xoff).astype('int')
for jj, (x_min, x_max) in enumerate(zip(xoff_bins_Q[:-1], xoff_bins_Q[1:])):
	sel = (Xoff >= x_min) & (Xoff < x_max) 
	coolness[sel] = jj+1

### the ones that have the most recent MM are most disturbed.
### we do 4 classes with equal number of each class
##time_of_last_MM = cosmo.age(1 / scale_of_last_MM - 1).value
##time_now = cosmo.age(zz).value
##delta_t_MM = time_now - time_of_last_MM
###
##delta_t_hist = n.histogram(delta_t_MM, bins=100)
##delta_t_fraction = n.hstack((0., n.cumsum(delta_t_hist[0]) / N_clu))
##delta_t_values_itp = interp1d(delta_t_hist[1], delta_t_fraction)
### if coolness is small, then it is disturbed
### if coolness is long, then it is relaxed
##coolness = delta_t_values_itp(delta_t_MM)


# implement correlated scatter for quantities
# as of now it is 100% correlated (false)

t = Table()

for col_name, unit_val in zip(f2[1].data.columns.names, f2[1].data.columns.units):
	t.add_column(Column(name=col_name, data=f2[1].data[col_name][cluster], unit=unit_val,dtype=n.float32 ) )

for col_name in f1[1].data.columns.names:
	if col_name == 'Mvir' or col_name == 'M200c' or col_name == 'M500c':
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster]/h, unit='', dtype=n.float32 ) )
	else:
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster], unit='', dtype=n.float32 ) )

t.add_column(Column(name='HALO_lineID', data=halo_lineID.astype('int'), unit='', dtype=n.int64 ) )
t.add_column(Column(name='CLUSTER_LX_soft_RF', data=n.log10(LX_out), unit='erg/s', dtype=n.float32))
t.add_column(Column(name='CLUSTER_kT', data=KT_OUT, unit='keV', dtype=n.float32))
t.add_column(Column(name='CLUSTER_LX_soft', data=n.log10(LX_obsF), unit='erg/s', dtype=n.float32))
t.add_column(Column(name='CLUSTER_FX_soft', data=CLU_FX_soft, unit='erg/(cm*cm*s)', dtype=n.float32))
t.add_column(Column(name='CLUSTER_coolness', data=coolness, unit='' , dtype=n.float32) )
t.add_column(Column(name='XRAY_image_path', data=n.array(path_2_images), unit='' , dtype='S22' ) )

for col_name, unit_val in zip(f3[1].data.columns.names, f3[1].data.columns.units):
	t.add_column(Column(name='galaxy_'+col_name, data=f3[1].data[col_name][cluster], unit='', dtype=n.float32 ) )

t.write(path_2_CLU_file, overwrite=True)

sys.exit()

# writes the results
print('writes results', time.time() - t0)
f = h5py.File(path_2_CLU_file, "a")
f.attrs['file_name'] = os.path.basename(path_2_CLU_file)
f.attrs['creator'] = 'JC'

# writes the results
halo_data = f.create_group('CLUSTERS')

halo_data.create_dataset('ids_cluster', data=ids_cluster)

ds = halo_data.create_dataset('LX_soft_RF', data=n.log10(LX_soft_RF))
ds.attrs['units'] = 'log10(L_X/[0.5-2keV, restFrame, erg/s])'

ds = halo_data.create_dataset('LX_soft', data=n.log10(LX_out))
ds.attrs['units'] = 'log10(L_X/[0.5-2keV, obsFrame, erg/s])'

halo_data.create_dataset('kT', data=KT_OUT)
ds.attrs['units'] = 'kT cin [keV]'

halo_data.create_dataset('FX_soft', data=CLU_FX_soft)
ds.attrs['units'] = 'F_X / [0.5-2keV, erg/cm2/s]'

halo_data.create_dataset('FX_soft_attenuated', data=CLU_FX_soft)
ds.attrs['units'] = 'F_X / [0.5-2keV, erg/cm2/s]'

halo_data.create_dataset('detected', data=detected)

halo_data.create_dataset('coolness', data=coolness)
ds.attrs['units'] = '0:very disturbed, 1:relaxed'

f.close()

### NOW THE FIGURES ###
if make_figure:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': 14})
    import matplotlib.pyplot as p

    fig_dir = os.path.join(
        os.environ['GIT_AGN_MOCK'],
        'figures',
        env,
        'clusters',
    )
    if os.path.isdir(fig_dir) == False:
        os.system('mkdir -p ' + fig_dir)

    fig_out = os.path.join(
        fig_dir,
        'Cluster_LX_soft_cin_vs_mass_' +
        baseName +
        '.png')
    X = M500c[cluster]
    Y = n.log10(LX_BB_18_cin) + scatter_1 * 0.27
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k+', rasterized=True)
    p.title(baseName)
    p.xscale('log')
    p.xlabel('M500c')
    p.ylabel('LX_soft cin')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(
        fig_dir,
        'Cluster_LX_soft_cex_vs_mass_' +
        baseName +
        '.png')
    X = M500c[cluster]
    Y = n.log10(LX_BB_18_cex) + scatter_1 * 0.27
    p.figure(1, (6., 5.5))
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
N_Mgas = 31.92  # Dietrich 18
N_kT = 2.18
N_L = 103.7
N_Lce = 102.66

slope_E_Mgas = 0.05  # Dietrich 18
slope_E_kT = 0.61
slope_E_L = 1.20
slope_E_Lce = 1.82

slope_M500_Mgas = 1.398  # Dietrich 18
slope_M500_kT = 0.66
slope_M500_L = 1.43  # 1.26*(1.+0.33*0.43)
slope_M500_Lce = 1.36  # 1.06*(1.+0.33*0.88)

scatter_Mgas = 0.106  # Dietrich 18
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

E035 = cosmo.efunc(0.35)

# converts logM500 to clusters observables


def m500_to_qty(logM500, z, slope_efunc, slope_m500, normalization): return n.e**normalization * \
    (cosmoMD.efunc(z) / E035)**(slope_efunc) * (10**(logM500 - n.log10(6) - 14))**(slope_m500)


def logM500_to_logMgas(
    logM500,
    z): return m500_to_qty(
        logM500,
        z,
        slope_E_Mgas,
        slope_M500_Mgas,
    N_Mgas)


def logM500_to_kT(
    logM500,
    z): return m500_to_qty(
        logM500,
        z,
        slope_E_kT,
        slope_M500_kT,
    N_kT)


def logM500_to_L(
    logM500,
    z): return m500_to_qty(
        logM500,
        z,
        slope_E_L,
        slope_M500_L,
    N_L)


def logM500_to_Lce(
    logM500,
    z): return m500_to_qty(
        logM500,
        z,
        slope_E_Lce,
        slope_M500_Lce,
    N_Lce)



#Mean_Mgas = n.log10(logM500_to_logMgas(log_m500c, z))
#V_scatter_Mgas = norm.rvs(loc=0, scale=scatter_Mgas, size=nCluster)
#VAL_Mgas = n.zeros_like(m500c_i)
#VAL_Mgas[ok] = Mean_Mgas + V_scatter_Mgas

#Mean_kT = logM500_to_kT(log_m500c, z)
#V_scatter_kT = norm.rvs(loc=0, scale=scatter_kT, size=nCluster)
#VAL_kT = n.zeros_like(m500c_i)
#VAL_kT[ok] = Mean_kT + V_scatter_kT

#Mean_L = n.log10(logM500_to_L(log_m500c, z))
#V_scatter_L = norm.rvs(loc=0, scale=scatter_L, size=nCluster)
#VAL_L = n.zeros_like(m500c_i)
#VAL_L[ok] = Mean_L + V_scatter_L

#Mean_Lce = n.log10(logM500_to_Lce(log_m500c, z))
#V_scatter_Lce = norm.rvs(loc=0, scale=scatter_Lce, size=nCluster)
#VAL_Lce = n.zeros_like(m500c_i)
#VAL_Lce[ok] = Mean_Lce + V_scatter_Lce
