"""
What it does
------------

Computes the X-ray cluster model in each shell

References
----------

comparat et al. in prep

Command to run
--------------

python3 004_0_cluster.py environmentVAR fileBasename image_boolean

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

image_boolean: yes or no will write or not the images for sixte / simput simulations

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, extinction, matplotlib


244 objects : 

182 XXL, 
49 HIFLUGS
12 XCOP
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
print('Creates the cluster file per shell')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


# import all pathes

env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "all_0.62840"
is_writing_images = sys.argv[3] # 'yes'
pixel_size_image = float(sys.argv[4]) # 20
z_snap = 1./float(baseName.split('_')[1])-1.
aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
print(env, baseName,z_snap)

# no need fo bias, it is already calibrated in the data
b_HS = float(sys.argv[5]) # 1.0 # 0.8
cov_mat_option = sys.argv[6] 

# initializes pathes to files
test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_CLU_file = os.path.join(test_dir, baseName + '_CLU_b'+str(int(b_HS*10))+'_CM_'+cov_mat_option + "_pixS_" + str(pixel_size_image)+'.fits')

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
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.) /omega(zz)
rho_critical = (cosmo.critical_density(z_snap).to(u.Msun*u.kpc**(-3))  ).value

f1 = fits.open(path_2_light_cone)
N_obj = len(f1[1].data['M500c'])
#cluster = ( n.log10(f1[1].data['M500c'] / h) > n.log10(1e13) ) # & ( f1[1].data['pid'] == -1 )
# 2e13 Msun
if env == "MD04":
	cluster = ( n.log10(f1[1].data['M500c'] / h) > 13.0 ) #& ( f1[1].data['pid'] == -1 )
	M500c = f1[1].data['M500c'][cluster] / h
elif env == "MD40" or env == "MD10" :
	cluster = ( n.log10(f1[1].data['M500c'] / h) > 13.7 ) #& ( f1[1].data['pid'] == -1 )
	M500c = f1[1].data['M500c'][cluster] / h
else:
	# in the UNIT simulation M500, M200 for sub haloes have had issues, so we replace them with an average rescaling from Mvir (Mvir for sub haloes looks ok)
	# this is a correction done for only <2% of the total number of haloes, which are less massive than their host and thus will most likely be outshine
	cluster_distinct = ( n.log10(f1[1].data['M500c'] / h) > 13.7 ) & ( f1[1].data['pid'] == -1 )
	cluster_sub = ( n.log10(0.6 * f1[1].data['Mvir']  / h) > 13.7 ) & ( f1[1].data['pid'] >= 0 )
	M500c = f1[1].data['M500c'][ cluster_distinct | cluster_sub ] / h
	M500c [cluster_sub] = 0.6 * f1[1].data['Mvir'][ cluster_sub ] / h
	cluster = (cluster_distinct) | (cluster_sub)

	
Mvir = f1[1].data['Mvir'][cluster] / h
Rvir = f1[1].data['Rvir'][cluster]
R500c = (DeltaVir_bn98(z_snap)/500. * M500c / Mvir)**(1./3.)*Rvir
# other formulae :
# R500c = n.power(3.* M500c / (4. * n.pi * 500. * rho_critical) ,1./3.)

logM500c = n.log10(M500c)
log_vmax = n.log10(f1[1].data['vmax'][cluster])
scale_of_last_MM = f1[1].data['scale_of_last_MM'][cluster]
N_clu = len(Mvir)
print(N_clu, 'clusters')
Xoff = f1[1].data['Xoff'][cluster] / f1[1].data['Rvir'][cluster]
NN,BB=n.histogram(Xoff, bins=1000)
NN_C = n.cumsum(NN)*1./N_clu
# from Xoff assign percentile
itp_xoff_cdf = interp1d( n.hstack((BB)), n.hstack((0.,NN_C)))
#from percentile assign Xoff
itp_xoff = interp1d( n.hstack((0.,NN_C)), n.hstack((BB)))
Xoff_PC = itp_xoff_cdf(Xoff)

b_to_a_500c = f1[1].data['b_to_a_500c'][cluster]
#f1.close()
halo_lineID = n.arange(N_obj)[cluster]
#ids_cluster = (halo_lineID*1000000000000 + abs(f1[1].data['id'][cluster]) ).astype('int64')
ids_cluster_str = n.array([aexp_str+str(el).zfill(16) for el in f1[1].data['id'][cluster]])
#aexp_str

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
#vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value)
#DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
#print('volume', vol, 'Mpc3', time.time() - t0)

f3 = fits.open(path_2_galaxy_file)
# mass = f3[1].data['SMHMR_mass'][cluster]  # log of the stellar mass
# f3.close()

path_2_cbp = os.path.join(os.environ['GIT_CBP'])

Mpc=3.0856776e+24
msun=1.98892e33
if N_clu < 100:
	nsim = N_clu*1000
else:
	nsim = N_clu*100

if cov_mat_option =='0':
	covor=n.loadtxt(os.path.join(path_2_cbp, 'covmat_xxl_hiflugcs_xcop.txt'))
	xgrid_ext=n.loadtxt(os.path.join(path_2_cbp, 'radial_binning.txt'))
	mean_log=n.loadtxt(os.path.join(path_2_cbp, 'mean_pars.txt'))
	coolfunc=n.loadtxt(os.path.join(path_2_cbp, 'coolfunc.dat'))
else :
	covor=n.loadtxt(os.path.join(path_2_cbp, 'covmat_'+cov_mat_option+'.txt'))
	xgrid_ext=n.loadtxt(os.path.join(path_2_cbp, 'radial_binning.txt'))
	mean_log=n.loadtxt(os.path.join(path_2_cbp, 'mean_pars_'+cov_mat_option+'.txt'))
	coolfunc=n.loadtxt(os.path.join(path_2_cbp, 'coolfunc.dat'))

def calc_lx(prof,kt,m5,z,fraction=1):
	"""
	Compute the X-ray luminosity in the profile
	to be extended to 3x r500c
	
	Compute r_{500c} :
		r_{500c} = \left(\frac{3 M_{500c}}{ 4. \pi 500 \rho_c(z)  }\right)^{1/3} [ cm ].
	
	profile_emission = profile x rescale_factor
	rescale_factor = \sqrt(kT/10.0) E^3(z)
	CF(kT) = cooling function, show the curve
	L_X(r) = \Sigma_{<r}( profile_emission r_{500c}^2 2 \pi x CF(kT) Mpc=3.0856776e+24 dx )
	L_{500c} = L_X(1)
	
	"""
	ez2 = cosmo.efunc(z)**2
	rhoc = cosmo.critical_density(z).value
	r500 = n.power(m5*msun/4.*3./n.pi/500./rhoc,1./3.)
	resfact = n.sqrt(kt/10.0)*n.power(ez2,3./2.)
	prof_em = prof * resfact # emission integral
	tlambda = n.interp(kt,coolfunc[:,0],coolfunc[:,1]) # cooling function
	dx = n.empty(len(xgrid_ext))
	dx[0]=xgrid_ext[0]
	dx[1:len(xgrid_ext)]=(n.roll(xgrid_ext,-1)-xgrid_ext)[:len(xgrid_ext)-1]
	#print(prof_em*xgrid_ext*r500**2*2.*n.pi*tlambda*Mpc*dx)
	lxcum = n.cumsum(prof_em*xgrid_ext*r500**2*2.*n.pi*tlambda*Mpc*dx) # riemann integral
	lx_500=n.interp(fraction,xgrid_ext,lxcum) # evaluated at R500
	return lx_500

profs=n.exp(n.random.multivariate_normal(mean_log,covor,size=nsim))

allz_i  = profs[:,len(mean_log)-3]
allkt_i = profs[:,len(mean_log)-1]
allm5_i = profs[:,len(mean_log)-2]

profiles_i = profs[:,:len(xgrid_ext)]

in_zbin = (allz_i<zmax+0.05)&(allz_i>zmin-0.05)

allz     = allz_i    [in_zbin]
allkt    = allkt_i   [in_zbin]
allm5    = allm5_i   [in_zbin]
profiles = profiles_i[in_zbin]
nsim2 = len(allz)

EM0 = -n.log10(profiles.T[0])
EM0_norm = ( EM0 - EM0.min() ) / n.max(EM0 - EM0.min())
NN,BB=n.histogram(EM0, bins=1000)
NN_C = n.cumsum(NN)*1./len(EM0)
# from EM0 assign percentile
itp_EM0_cdf = interp1d( n.hstack((BB)), n.hstack((0.,NN_C)))
#from percentile assign Xoff
itp_EM0 = interp1d( n.hstack((0.,NN_C)), n.hstack((BB)))
EM0_PC = itp_EM0_cdf(EM0)

alllx           = n.empty(nsim2)
alllx_halfr500  = n.empty(nsim2)
alllx_twicer500 = n.empty(nsim2)
for i in range(nsim2):
	tprof=profiles[i]
	alllx[i]          =calc_lx(tprof,allkt[i],allm5[i],allz[i], 1)
	alllx_halfr500[i] =calc_lx(tprof,allkt[i],allm5[i],allz[i], 0.5)
	alllx_twicer500[i]=calc_lx(tprof,allkt[i],allm5[i],allz[i], 2.0)

print(len(allz),'profiles created', time.time()-t0)
from sklearn.neighbors import BallTree

Tree_profiles = BallTree(n.transpose([allz, n.log10(allm5), EM0_PC]))
#DATA = n.array([[z_i, m_i] for z_i, m_i in zip(zzs, logM500c) ])
DATA = n.transpose([zzs, n.log10(M500c*b_HS), Xoff_PC])
ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)

ids = n.hstack((ids_out))

# LX_out 0.5-2 rest frame erg/s
# to convert to 0.5-2 observed frame 
LX_out = alllx[ids]
LX_out_halfr500  = alllx_halfr500 [ids] 
LX_out_twicer500 = alllx_twicer500[ids]
KT_OUT = allkt[ids]
CBP_redshift = allz[ids]
CBP_M500c    = n.log10(allm5)[ids]
CBP_EM0      = itp_EM0(EM0_PC[ids])
print('profiles matched', time.time()-t0)

# maximum kT for the k-correction
KT_OUT[KT_OUT>19.4] = 19.4

# attenuation grid should cover down to 0.05 keV 
itp_logNH, itp_z, itp_kt, itp_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters.txt"), unpack=True )

nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
#z_vals = n.hstack(( n.arange(0.,0.7,0.05), n.arange(0.8, 4.5, 0.1)))#[0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6] ))
#kT_vals = n.hstack(( n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))
#kT_vals = n.hstack(([0.1, 0.2], n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))
kT_vals = 10**n.arange(-1,1.3,0.01)
z_vals = n.hstack((n.array([0., 0.01]), 10**n.arange(n.log10(0.02), n.log10(4.), 0.01)))

XX_nh, YY_z, ZZ_kt = n.meshgrid(nh_vals, z_vals, kT_vals)

shape_i = XX_nh.shape

matrix_z_nh_kt = itp_frac_obs.reshape(shape_i)

from scipy.interpolate import RegularGridInterpolator
attenuation_3d = RegularGridInterpolator((z_vals, n.log10(nh_vals*1e22), kT_vals), matrix_z_nh_kt)

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
pix_ids = healpy.ang2pix( 512,    n.pi /    2. -    g_lat *    n.pi /    180.,    g_lon *    n.pi /    180.,    nest=True)
flux_lim_data = fits.open(path_2_flux_limits)
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3
flux_limit = flux_lim_data[1].data['flux_limit_eRASS8'][pix_ids]
# HPX 8
#HEALPIX_8 = healpy.ang2pix(8, n.pi/2. - dec*n.pi/180. , ra*n.pi/180. , nest=True)

print('flux limit file opened', time.time() - t0)
#if is_writing_images=="yes" and pixel_size_image>=18:
	#detected = (CLU_FX_soft > 10**(flux_limit)) 
#elif is_writing_images=="yes" and pixel_size_image<18 :
	#detected = (CLU_FX_soft > 10**(-15)) & (HEALPIX_8==0) 
#else :
	#detected = (CLU_FX_soft > 10**(-15)) # create way too many points
	
# relaxation state
# based on MD40, quartiles of Xoff
#xoff_bins_Q = n.array([0.,  4.28608992e-02, 6.59490117e-02, 1.00837135e-01, 5.7516e-01])
#coolness = n.ones_like(Xoff).astype('int')
#for jj, (x_min, x_max) in enumerate(zip(xoff_bins_Q[:-1], xoff_bins_Q[1:])):
	#sel = (Xoff >= x_min) & (Xoff < x_max) 
	#coolness[sel] = jj+1

dir_2_SMPT_image = os.path.join(os.environ[env], "cat_CLU_SIMPUT_b" + str(int(b_HS*10)) + '_CM_' + cov_mat_option + "_pixS_" + str(pixel_size_image), 'cluster_images', baseName )
if os.path.isdir(dir_2_SMPT_image) == False:
    os.system('mkdir -p ' + dir_2_SMPT_image)

print('creates SIMPUT images')
path_2_images = []
# writes images using these profiles

# puts aminimum to ellipticity
# otherwise images cannot be computed
b_to_a_500c[b_to_a_500c<0.01] = 0.01

angularSize_per_pixel = n.zeros_like(ra)
#for jj, (profile_i, r500c_i, conversion_arcmin, file_name, b_a, det) in enumerate(zip(profiles[ids], R500c,  cosmo.kpc_proper_per_arcmin(zz).value, ids_cluster_str, b_to_a_500c, detected)):
for jj, (profile_i, r500c_i, conversion_arcmin, file_name, b_a) in enumerate(zip(profiles[ids], R500c,  cosmo.kpc_proper_per_arcmin(zz).value, ids_cluster_str, b_to_a_500c)):
	image_file = os.path.join(dir_2_SMPT_image, file_name+'.fits')
	#print(image_file)
	path_2_images.append(os.path.join('cluster_images', baseName, file_name) )
	if is_writing_images == 'yes' :
		#if det :
		x_coord = n.hstack(( 0., xgrid_ext*r500c_i/conversion_arcmin, 10*xgrid_ext[-1]*r500c_i/conversion_arcmin, 10000000 )) 
		y_coord = n.hstack(( profile_i[0], profile_i, 0., 0. ))
		profile = interp1d(x_coord, y_coord)
		truncation_radius = 2 * r500c_i/ conversion_arcmin
		n_pixel = 2*int(truncation_radius*60/pixel_size_image)+1
		matrix, angularSize_per_pixel_j = create_matrix(profile, n_pixel = n_pixel, b_a = b_a, truncation_radius = truncation_radius)
		angularSize_per_pixel[jj] = angularSize_per_pixel_j
		print(r500c_i, conversion_arcmin, r500c_i/ conversion_arcmin, 60*angularSize_per_pixel_j, n_pixel, image_file)
		if os.path.isfile(image_file)==False:
			write_img(matrix, image_file, n_pixel = n_pixel, angularSize_per_pixel=angularSize_per_pixel_j)

print('angularSize_per_pixel', angularSize_per_pixel)

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
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster]/h, unit='Msun', description='halo mass in aperture. No little h', dtype=n.float32 ) )
	elif col_name == 'id' or col_name == 'pid' :
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster], unit='', dtype=n.int64 ) )
	else:
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster], unit='', dtype=n.float32 ) )

t.add_column(Column(name='HALO_lineID', data=halo_lineID.astype('int'), unit='', dtype=n.int64 ) )
t.add_column(Column(name='CLUSTER_LX_soft_RF', data=n.log10(LX_out), unit='', description='rest-frame log10(L_X [erg/s]), R < R500c, 0.5-2keV', dtype=n.float32))
t.add_column(Column(name='CLUSTER_LX_soft_RF_halfR500', data=n.log10(LX_out_halfr500), unit='', description='rest-frame log10(L_X [erg/s]), R < 0.5 R500c, 0.5-2keV', dtype=n.float32))
t.add_column(Column(name='CLUSTER_LX_soft_RF_twiceR500', data=n.log10(LX_out_twicer500), unit='', description='rest-frame log10(L_X [erg/s]), R < 2 R500c, 0.5-2keV', dtype=n.float32))
t.add_column(Column(name='CLUSTER_kT', data=KT_OUT, unit='keV', description='Cluster temperature', dtype=n.float32))
t.add_column(Column(name='CLUSTER_LX_soft', data=n.log10(LX_obsF), unit='', description='Observed cluster luminosity log10(L_X [erg/s]), R < R500c, 0.5-2keV', dtype=n.float32))
t.add_column(Column(name='CLUSTER_FX_soft', data=CLU_FX_soft, unit='erg/(cm*cm*s)', description='Observed cluster flux log10(F_X [erg/s]), R < R500c, 0.5-2keV', dtype=n.float32))
t.add_column(Column(name='XRAY_image_path', data=n.array(path_2_images), unit='', description='relative path to the image for sixte simulation' , dtype='S60' ) )
t.add_column(Column(name='angularSize_per_pixel', data=angularSize_per_pixel, unit='arcmin' , description='angular size in arcminute of a pixel in the image XRAY_image_path', dtype=n.float32 ) )
t.add_column(Column(name='CBP_redshift', data=CBP_redshift, unit='', description='redshift of the matched simulated profile' , dtype=n.float32 ) )
t.add_column(Column(name='CBP_M500c', data=CBP_M500c, unit='', description='M500c of the matched simulated profile' , dtype=n.float32 ) )
t.add_column(Column(name='CBP_EM0', data=CBP_EM0, unit='', description='EM(0) of the matched simulated profile' , dtype=n.float32 ) )
t.add_column(Column(name='flux_limit_eRASS8_pt', data=flux_limit, unit='', description='Realization of the point source flux limit from Clerc et al. 2018 log10(F_X [erg/s]) 0.5-2keV' , dtype=n.float32 ) )

for col_name, unit_val in zip(f3[1].data.columns.names, f3[1].data.columns.units):
	t.add_column(Column(name='galaxy_'+col_name, data=f3[1].data[col_name][cluster], unit='', dtype=n.float32 ) )

t.write(path_2_CLU_file, overwrite=True)
