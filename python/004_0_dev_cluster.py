"""
What it does
------------

Adds the X-ray cluster model on a light cone

Randomly samples clusters profiles and temperatures from the recent method by dev. D. Eckert :
https://github.com/ygolomsochtiwsretsulc/cluster-brightness-profiles

Assigns cluster properties to the closest in the (z, M500c) plane

References
----------

Eckert et al. in prep

Uses data from XXL, HIFLUGS, XCOP, 244 objects : 
  * 182 XXL 
  * 49 HIFLUGS
  * 12 XCOP

Upcoming improvement, add ellipticities

Command to run
--------------

python3 004_0_dev_cluster.py 

arguments
---------

should add in argument 
 * path to the input light cone file
 * path to the output light cone file

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, extinction, matplotlib, xspec

https://github.com/ygolomsochtiwsretsulc/cluster-brightness-profiles

H1 maps 

To do 
-----



"""
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
from astropy.table import Table, Column
from scipy.linalg import eigh, cholesky
from scipy.stats import norm

import h5py
import numpy as n
print('Creates the h5 Cluster file per shell')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

# import all pathes
env = 'MD10' 
make_figure = True

# initializes pathes to files
path_2_light_cone = os.path.join(os.environ[env], 'MD10_eRO_CLU.fit')
path_2_result = os.path.join(os.environ[env], 'MD10_eRO_CLU_Eckert19.fit')
data = Table.read(path_2_light_cone)#[1].data

data.remove_columns(['LX_soft_cin','LX_soft_cex','TX_cin','TX_cex','MICM_cex','FX_soft','FX_soft_attenuated','detected','scatter_1','coolness'])

## x ray extinction for clusters, from A. Finoguenov, due to foreground
#path_2_NH_law = os.path.join(
	#os.environ['GIT_AGN_MOCK'],
	#"data",
	#"xray_k_correction",
	#'nh_flux.tbl')
#NH_DATA = n.loadtxt(path_2_NH_law, unpack=True)
#nh_law = interp1d(
	#n.hstack(
		#(-1.e30,
		#1.0 * 10**10,
		#NH_DATA[0],
		#1.0 * 10**30)),
	#n.hstack(
			#(1.,
			#1.,
			#NH_DATA[1][0] / NH_DATA[1],
			#7.03612982e-09)))

# eRosita flux limits for point sources
path_2_flux_limits = os.path.join(
	os.environ['GIT_AGN_MOCK'],
	"data",
	"erosita",
	"flux_limits.fits")

# simulation setup
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


area = 27143./2. # deg2
sel_area = (abs(data['g_lat'])>20) # &(data['dec']<0)

Mvir = data['HALO_Mvir']
M500c = data['HALO_M500c'].data
logM500c = n.log10(M500c )
log_vmax = n.log10(data['HALO_vmax'])
b_2_a = data['HALO_b_to_a_500c']

HALO_rs = data['HALO_rs'] 
HALO_rvir = data['HALO_rvir'] 
zzs = data['redshift_S']

N_clu = len(Mvir)


path_2_cbp = os.path.join(os.environ['GIT_CBP'])

Mpc=3.0856776e+24
msun=1.98892e33
nsim = 1000000
covor=n.loadtxt(os.path.join(path_2_cbp, 'covmat_xxl_hiflugcs_xcop.txt'))
xgrid_ext=n.loadtxt(os.path.join(path_2_cbp, 'radial_binning.txt'))
mean_log=n.loadtxt(os.path.join(path_2_cbp, 'mean_pars.txt'))
coolfunc=n.loadtxt(os.path.join(path_2_cbp, 'coolfunc.dat'))

def calc_lx(prof,kt,m5,z):
	"""
	Compute the X-ray luminosity in the profile

	to be exteneded to 3x r500c
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



from sklearn.neighbors import BallTree

Tree_profiles = BallTree(n.transpose([allz, n.log10(allm5)]))
DATA = n.transpose([zzs, logM500c])

ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)

ids = n.hstack((ids_out))

# LX_out 0.5-2 rest frame erg/s
# to convert to 0.5-2 observed frame 
LX_out = alllx[ids]
KT_OUT = allkt[ids]

data.add_column(Column(name='LX_soft_RF', data=LX_out, unit='erg/s'))
data.add_column(Column(name='kT', data=KT_OUT, unit='keV'))


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

#fig_out = os.path.join(fig_dir, 'Cluster_LX_LX.png')

#xs = n.arange(42, 45, 0.1)
#p.figure(1, (6., 5.5))
##p.tight_layout()
#p.plot(n.log10(LX_out), data['LX_soft_cin'],' k,', rasterized=True)
#p.plot(xs, xs, 'r--')
##p.xscale('log')
##p.yscale('log')
#p.xlabel('LX Eckert19')
#p.ylabel('LX Mantz16')
#p.grid()
#p.savefig(fig_out)
#p.clf()


#fig_out = os.path.join(fig_dir, 'Cluster_kT_kT.png')

#xs = n.arange(0.5, 20, 0.1)
#p.figure(1, (6., 5.5))
##p.tight_layout()
#p.plot(KT_out, 10**data['TX_cin'],' k,', rasterized=True)
#p.plot(xs, xs, 'r--')
#p.xscale('log')
#p.yscale('log')
#p.xlabel('kT Eckert19')
#p.ylabel('kT Mantz16')
#p.grid()
#p.savefig(fig_out)
#p.clf()

# Compute K correction
# add X spectra with NH
# write 1 million images with ellipticity of the simulation

# write simput catalogues

##from scipy.interpolate import RegularGridInterpolator
##def f(x,y,z):
    ##return 2 * x**3 + 3 * y**2 - z
##x = np.linspace(1, 4, 11)
##y = np.linspace(4, 7, 22)
##z = np.linspace(7, 9, 33)
##data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))

##my_interpolating_function = RegularGridInterpolator((x, y, z), data)
##pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
##my_interpolating_function(pts)

itp_logNH, itp_z, itp_kt, itp_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters.txt"), unpack=True )

nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
z_vals = n.hstack(( n.arange(0.,0.7,0.05), [0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6] ))
kT_vals = n.hstack(( n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))

XX_nh, YY_z, ZZ_kt = n.meshgrid(nh_vals, z_vals, kT_vals)

shape_i = XX_nh.shape

matrix = itp_frac_obs.reshape(shape_i)

from scipy.interpolate import RegularGridInterpolator
attenuation_3d = RegularGridInterpolator((z_vals, n.log10(nh_vals*1e22), kT_vals), matrix)

sel = (data['galactic_NH']>1e19)
bad_NH = (sel==False)
data['galactic_NH'][bad_NH] = 1.000001e19

k_correction_3d = attenuation_3d( n.transpose([zzs, n.ones_like(n.log10(data['galactic_NH']))*20., KT_OUT]))

itp_z, itp_kt, itp_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters_no_nH.txt"), unpack=True )


z_vals = n.hstack(( n.arange(0.,0.7,0.05), [0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6] ))
kT_vals = n.hstack(( n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))

YY_z, ZZ_kt = n.meshgrid(z_vals, kT_vals)

shape_i = YY_z.shape

matrix_2d = itp_frac_obs.reshape(shape_i)

from scipy.interpolate import RegularGridInterpolator
attenuation_2d = RegularGridInterpolator((kT_vals, z_vals), matrix_2d)

k_correction_2d = attenuation_2d( n.transpose([KT_OUT, zzs]))


attenuate_X_logNH, attenuate_Y_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "nh_attenuation_clusters.txt"), unpack=True )

#att = attenuation( n.transpose([zzs, n.log10(data['galactic_NH']), KT_OUT]))

LX_obsF = LX_out / k_correction_2d 
dL_cm = data['dL_cm']
CLU_FX_soft_1 = LX_obsF / (4 * n.pi * dL_cm.data**2.)

itp_attenuation = interp1d(attenuate_X_logNH, attenuate_Y_frac_obs)
att = itp_attenuation(n.log10(data['galactic_NH']))
CLU_FX_soft = CLU_FX_soft_1 * att


data.add_column(Column(name='LX_soft', data=LX_obsF, unit='erg/s'))
data.add_column(Column(name='FX_soft', data=CLU_FX_soft, unit='erg/cm2/s'))

data.write(path_2_result)


sys.exit()


"""


topcat -stilts plot2plane \
   xpix=487 \
   ylog=true \
    xlabel='redshift_R / real space' \
    ylabel=LX_soft/LX_soft_RF grid=true \
    fontsize=16 fontstyle=serif \
    fontweight=bold \
   xmin=0 xmax=1.54 ymin=0.589 ymax=1.658 \
   auxmap=accent auxmin=0.98 auxmax=10 \
   auxvisible=true auxlabel='kT / keV' \
   legend=false \
   layer=Mark \
      in=/home/comparat/data/MultiDark/MD_1.0Gpc/MD10_eRO_CLU_Eckert19.fit \
      x=redshift_R y=LX_soft/LX_soft_RF \
       aux=kT \
      shading=aux \
   omode=out out=/home/comparat/software/lss_mock_dev/figures/MD10/clusters/k_correction.png 


"""







coolness = data['coolness']
ra =  data['ra']
dec = data['dec']
zz = data['redshift_R']
dL_cm = data['dL_cm']
galactic_NH = data['galactic_NH']
galactic_ebv = data['galactic_ebv']
g_lat = data['g_lat']
g_lon = data['g_lon']
ecl_lat = data['ecl_lat']
ecl_lon = data['ecl_lon']

#print('coordinate file opened', time.time() - t0)

# correlated random variables 
# https://scipy-cookbook.readthedocs.io/items/CorrelatedRandomSamples.html
# Values from Farahi et al. 2019 L_X, T_X
r = n.array([
        [  1, 0.49],
        [  0.49,  1],
    ])
# Generate samples from three independent normally distributed random
# variables (with mean 0 and std. dev. 1).
x = norm.rvs(size=(2, N_clu))
# We need a matrix `c` for which `c*c^T = r`.  We can use, for example,
# the Cholesky decomposition, or the we can construct `c` from the
# eigenvectors and eigenvalues.
# Compute the Cholesky decomposition.
c = cholesky(r, lower=True)
#evals, evecs = eigh(r)
# Construct c, so c*c^T = r.
#c = n.dot(evecs, n.diag(n.sqrt(evals)))
# Convert the data to correlated random variables. 
y = n.dot(c, x)
scatter_1, scatter_2 =y

# Mantz et al. 2016 scaling relations 

N_Mgas = 31.98 
N_kT 	= 2.18
N_L 	= 103.7
N_Lce 	= 102.66

slope_E_Mgas 	= -0.11 
slope_E_kT 	= 0.61
slope_E_L 	= 1.20
slope_E_Lce 	= 1.82

slope_M500_Mgas= 1.04
slope_M500_kT 	= 0.66
slope_M500_L 	= 1.26
slope_M500_Lce = 1.06

scatter_Mgas = 0.086
scatter_kT = 0.18
scatter_L = 0.24
scatter_Lce = 0.17

E035 = cosmo.efunc(0.35)

# converts logM500 to clusters observables
m500_to_qty = lambda logM500, z, slope_efunc, slope_m500, normalization : n.e**normalization * (cosmoMD.efunc(z)/E035)**(slope_efunc) * (10**(logM500-n.log10(6)-14))**(slope_m500)

#logM500_to_logMgas 	= lambda logM500, z : m500_to_qty( logM500, z, slope_E_Mgas, slope_M500_Mgas, N_Mgas)
logM500_to_kT 		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_kT, slope_M500_kT, N_kT)
logM500_to_L 		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_L, slope_M500_L, N_L)
#logM500_to_Lce		= lambda logM500, z : m500_to_qty( logM500, z, slope_E_Lce, slope_M500_Lce, N_Lce)

Mean_kT = logM500_to_kT(logM500c, zz)
VAL_kT = Mean_kT + scatter_1 * scatter_kT

Mean_L = n.log10(logM500_to_L(logM500c, zz))
VAL_L = Mean_L + scatter_2 * scatter_L

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

# Flux observed and attenuated
attenuation = nh_law(galactic_NH)
CLU_FX_soft = 10**VAL_L / (4 * n.pi * dL_cm**2.)
CLU_FX_soft_attenuated = attenuation * CLU_FX_soft
detected = (CLU_FX_soft_attenuated > 10**(flux_limit))

data['LX_soft_cin'] = VAL_L
data['TX_cin'] = VAL_kT

data.remove_column('LX_soft_cex')
data.remove_column('TX_cex')
data.remove_column('MICM_cex')
data.remove_column('scatter_1')

data.add_column(Column(name='scatter_kT', data=scatter_1, unit=''))
data.add_column(Column(name='scatter_LX', data=scatter_2, unit=''))

data.write(path_2_result)

sys.exit()

# writes the results
print('writes results', time.time() - t0)
f = h5py.File(path_2_CLU_file, "a")
f.attrs['file_name'] = os.path.basename(path_2_CLU_file)
f.attrs['creator'] = 'JC'

# writes the results
halo_data = f.create_group('CLUSTERS')

halo_data.create_dataset('ids_cluster', data=ids_cluster)

ds = halo_data.create_dataset(
    'LX_soft_cin',
    data=n.log10(LX_BB_18_cin) +
    scatter_1 *
    0.27)
ds.attrs['units'] = 'log10(L_X/[0.5-2keV, cin, erg/s])'

halo_data.create_dataset(
    'LX_soft_cex',
    data=n.log10(LX_BB_18_cex) +
    scatter_1 *
    0.27)
ds.attrs['units'] = 'log10(L_X/[0.5-2keV, cex, erg/s])'

halo_data.create_dataset(
    'TX_cin',
    data=n.log10(TX_BB_18_cin) +
    scatter_1 *
    0.179)
ds.attrs['units'] = 'log10(T cin [keV])'

halo_data.create_dataset(
    'TX_cex',
    data=n.log10(TX_BB_18_cex) +
    scatter_1 *
    0.128)
ds.attrs['units'] = 'log10(T cex [keV])'

halo_data.create_dataset(
    'MICM_cex',
    data=n.log10(MICM_BB_18_cex) +
    scatter_1 *
    0.098)
ds.attrs['units'] = 'log10(M_ICM [Msun]))'

halo_data.create_dataset('FX_soft', data=CLU_FX_soft)
ds.attrs['units'] = 'F_X / [0.5-2keV, erg/cm2/s]'

halo_data.create_dataset('FX_soft_attenuated', data=CLU_FX_soft_attenuated)
ds.attrs['units'] = 'F_X / [0.5-2keV, erg/cm2/s]'

halo_data.create_dataset('detected', data=detected)

halo_data.create_dataset('scatter_1', data=scatter_1)

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

E035 = cosmoMD.efunc(0.35)

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


z = f1.attrs['redshift']
m500c_i = f1['/halo_properties/M500c'].value
ok = (m500c_i > 1e13)
log_m500c = n.log10(m500c_i[ok])

nCluster = len(log_m500c)
#rds = (n.random.rand(len(log_m500c))-0.5)*2.

Mean_Mgas = n.log10(logM500_to_logMgas(log_m500c, z))
V_scatter_Mgas = norm.rvs(loc=0, scale=scatter_Mgas, size=nCluster)
VAL_Mgas = n.zeros_like(m500c_i)
VAL_Mgas[ok] = Mean_Mgas + V_scatter_Mgas

Mean_kT = logM500_to_kT(log_m500c, z)
V_scatter_kT = norm.rvs(loc=0, scale=scatter_kT, size=nCluster)
VAL_kT = n.zeros_like(m500c_i)
VAL_kT[ok] = Mean_kT + V_scatter_kT

Mean_L = n.log10(logM500_to_L(log_m500c, z))
V_scatter_L = norm.rvs(loc=0, scale=scatter_L, size=nCluster)
VAL_L = n.zeros_like(m500c_i)
VAL_L[ok] = Mean_L + V_scatter_L

Mean_Lce = n.log10(logM500_to_Lce(log_m500c, z))
V_scatter_Lce = norm.rvs(loc=0, scale=scatter_Lce, size=nCluster)
VAL_Lce = n.zeros_like(m500c_i)
VAL_Lce[ok] = Mean_Lce + V_scatter_Lce
