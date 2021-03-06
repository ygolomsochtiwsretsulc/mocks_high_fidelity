"""
What it does
------------

Create a file with galaxy properties.

Computes
 - stellar mass following Moster et al. 2013, 2018
 - star formation rates following the main sequence of galaxies following Whitaker et al 2012
 - star formation rate for quiescent galaxies following fits on COSMOS data, see Comparat et al. in prep
 - X-ray luminosities of galaxies after Aird et al. 2018
 - r-band magnitude in the SDSS AB system doing abundance matching with the luminosity function of Loveday et al. 2016

References
----------

Whitaker 2012, 2014
https://ui.adsabs.harvard.edu/abs/2012ApJ...754L..29W
https://ui.adsabs.harvard.edu/abs/2014ApJ...795..104W
Loveday et al. 2016

Ilber et al. 2013


Command to run
--------------

python3 002_0_galaxy.py environmentVAR fileBasename

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

It will find the input files :
$environmentVAR/hlists/fits/${fileBasename}.fits
$environmentVAR/hlists/fits/${fileBasename}_coordinates.fits

And will write outputs in h5 format:
$environmentVAR/hlists/fits/${fileBasename}_galaxy.fits

Figures (Optional)
if the variable 'make_figure' is set to True, then figures will be created in the git repo here :
$GIT_AGN_MOCK/figures/$environmentVAR/galaxy/
It makes relevant histograms and scatter plots for all columns present in the new file.

cd /data39s/simulation_2/MD/MD_1.0Gpc/UniverseMachine
rsync -avz joco@draco01.mpcdf.mpg.de:~/ptmp_joco/simulations/MD_1.0Gpc/UniverseMachine/sfr_catalog_0.7* .

cd /data39s/simulation_2/MD/MD_1.0Gpc/UniverseMachine
rsync -avz joco@draco01.mpcdf.mpg.de:~/ptmp_joco/simulations/MD_1.0Gpc/UniverseMachine/sfr_catalog_0.8* .

cd /data39s/simulation_2/MD/MD_1.0Gpc/UniverseMachine
rsync -avz joco@draco01.mpcdf.mpg.de:~/ptmp_joco/simulations/MD_1.0Gpc/UniverseMachine/sfr_catalog_0.9* .

cd /data39s/simulation_2/MD/MD_1.0Gpc/UniverseMachine
rsync -avz joco@draco01.mpcdf.mpg.de:~/ptmp_joco/simulations/MD_1.0Gpc/UniverseMachine/sfr_catalog_1* .

wget http://behroozi.users.hpc.arizona.edu/UniverseMachine/DR1/SMDPL_SFR/sfr_catalog_0.700300.bin .

cd $MD04/UniverseMachine
rsync -avz joco@draco01.mpcdf.mpg.de:~/ptmp_joco/simulations/MD_0.4Gpc/UniverseMachine/sfr_catalog_* .

"""
import sys, os, scipy
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp1d
import h5py
import astropy.io.fits as fits
from astropy.table import Table, Column
import numpy as n
import time
print('Adds galaxy properties ')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

dtype_UM = n.dtype(
	dtype=
	[('id', 'i8'),('descid','i8'),('upid','i8'),
	('flags', 'i4'), ('uparent_dist', 'f4'),
	('pos', 'f4', (6)), ('vmp', 'f4'), ('lvmp', 'f4'),
	('mp', 'f4'), ('m', 'f4'), ('v', 'f4'), ('r', 'f4'),
	('rank1', 'f4'), ('rank2', 'f4'), ('ra', 'f4'),
	('rarank', 'f4'), ('A_UV', 'f4'), ('sm', 'f4'), 
	('icl', 'f4'), ('sfr', 'f4'), ('obs_sm', 'f4'), 
	('obs_sfr', 'f4'), ('obs_uv', 'f4'), ('empty', 'f4')],
	align=True)

# ID: Unique halo ID
# DescID: ID of descendant halo (or -1 at z=0).
# UPID: -1 for central halos, otherwise, ID of largest parent halo
# Flags: Ignore
# Uparent_Dist: Ignore
# pos[6]: (X,Y,Z,VX,VY,VZ)
# X Y Z: halo position (comoving Mpc/h)
# VX VY VZ: halo velocity (physical peculiar km/s)
# M: Halo mass (Bryan & Norman 1998 virial mass, Msun/h)
# V: Halo vmax (physical km/s)
# MP: Halo peak historical mass (BN98 vir, Msun/h)
# VMP: Halo vmax at the time when peak mass was reached.
# R: Halo radius (BN98 vir, comoving kpc/h)
# Rank1: halo rank in Delta_vmax (see UniverseMachine paper)
# Rank2, RA, RARank: Ignore
# A_UV: UV attenuation (mag)
# SM: True stellar mass (Msun)
# ICL: True intracluster stellar mass (Msun)
# SFR: True star formation rate (Msun/yr)
# Obs_SM: observed stellar mass, including random & systematic errors (Msun)
# Obs_SFR: observed SFR, including random & systematic errors (Msun/yr)
# Obs_UV: Observed UV Magnitude (M_1500 AB)


env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "sat_0.62840"
aexp_str = baseName.split('_')[1]
print(env, baseName, aexp_str)
make_figure = True
make_figure = False


def get_a(baseName):
    alp = baseName.split('_')[1]
    print('a=', alp)
    return float(alp)

a_snap = get_a(baseName)

# import all pathes
test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_universeMachine_file = os.path.join(os.environ[env], 'UniverseMachine', 'sfr_catalog_' + aexp_str + '0.bin' )

if a_snap < 0.40:
    fraction = 0.3
else:
    fraction = 1.0

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

print('opens light cone ')
f1 = fits.open(path_2_light_cone)
f1_id = f1[1].data['id']
Mvir = f1[1].data['Mvir'] / h
log_vmax = n.log10(f1[1].data['vmax'])
N_obj = len(Mvir)
f1.close()

print('opens coordinates ')
f2 = fits.open(path_2_coordinate_file) # [1].data
zz = f2[1].data['redshift_R']
N_halos = len(zz)
f2.close()

print('opens Universe Machine ')
um_galaxies = n.fromfile(path_2_universeMachine_file, dtype = dtype_UM)

# A_UV    : UV attenuation (mag)
# SM      : True stellar mass (Msun)
# ICL     : True intracluster stellar mass (Msun)
# SFR     : True star formation rate (Msun/yr)
# Obs_SM  : observed stellar mass, including random & systematic errors (Msun)
# Obs_SFR : observed SFR, including random & systematic errors (Msun/yr)
# Obs_UV  : Observed UV Magnitude (M_1500 AB)

snap_id = um_galaxies['id']
A_UV    = um_galaxies['A_UV']    
SM      = um_galaxies['sm']      
ICL     = um_galaxies['icl']     
SFR     = um_galaxies['sfr']     
Obs_SM  = um_galaxies['obs_sm']  
Obs_SFR = um_galaxies['obs_sfr'] 
Obs_UV  = um_galaxies['obs_uv']  

print('matches IDs')
snap_id_sortIDX = n.argsort(snap_id)

# link DM halo ids and UM ids
uniq_ids = n.unique(f1_id)#.astype('int32')) 

out = n.searchsorted(snap_id[snap_id_sortIDX], uniq_ids)
# 1st verification: must be true everywhere
print( len(snap_id[snap_id_sortIDX[out]]), len((snap_id[snap_id_sortIDX[out]]==uniq_ids).nonzero()[0]), snap_id[snap_id_sortIDX[out]]==uniq_ids )
# extract the indexes ?

id_map = {}
for ii,jj in zip(uniq_ids, snap_id_sortIDX[out]):
	id_map[ii]=jj

id_2_UM = n.array([id_map[el] for el in f1_id])

UM_snap_id = snap_id[id_2_UM]
# 2nd verification: must be true everywhere
print( len(f1_id), len((f1_id==UM_snap_id).nonzero()[0]), f1_id==UM_snap_id )
UM_A_UV    = A_UV   [id_2_UM]
UM_SM      = SM     [id_2_UM]
UM_ICL     = ICL    [id_2_UM]
UM_SFR     = SFR    [id_2_UM]
UM_Obs_UV  = Obs_UV [id_2_UM]
mass       = Obs_SM [id_2_UM]
sfr        = Obs_SFR[id_2_UM]

#N_galaxies = len(zz)
#volume = (cosmo.comoving_volume(n.max(zz)) - cosmo.comoving_volume(n.min(zz))).value
#print(volume)

# STELLAR MASS
# Equations 1 of Comparat et al. 2019

# def meanSM(Mh, z): return n.log10(Mh * 2. * (0.0351 - 0.0247 * z / (1. + z)) / ((Mh / (10**(11.79 + 1.5 * z / (1. + z))))** (- 0.9 + 0.5 * z / (1. + z)) + (Mh / (10**(11.79 + 1.5 * z / (1. + z))))**(0.67 + 0.2 * z / (1. + z))))

#mean_SM = meanSM(Mvir, zz)
#rds = norm.rvs(loc=0, scale=0.15, size=N_halos)
#mass = mean_SM + rds  # fun(mean_SM)
print('masses', mass[:10], time.time() - t0)

# quiescent fraction fitted on COSMOS, Ilbert et al. 2013 v2 catalogue
def scatter_Qf(z): return - 0.45 * (1 + z) + 1.54


def log10_M0_Qf(z): return 9.71 + 0.78 * (1 + z)


def fraction_Qf(mass, z): return 0.5 + 0.5 * \
    erf((mass - log10_M0_Qf(z)) / scatter_Qf(z))


rds_Qf = n.random.random(N_obj)

frac = fraction_Qf(n.log10(mass), zz)
frac[ zz > 2 ] = 0.
SF = (rds_Qf > frac)
QU = (SF == False)

# Hard X-ray emission, after Aird et al. 2018
def galaxy_lx(redshift, mass, sfr):
	return 10**(28.81) * (1 + redshift)**(3.9) * mass + \
        10**(39.5) * (1 + redshift)**(0.67) * sfr**(0.86)


gal_LX = galaxy_lx(zz, mass, sfr)
print('gal_LX', gal_LX[:10])

t = Table()
t.add_column(Column(name='SMHMR_mass', data=n.log10(mass), unit='log10(Msun)'))
t.add_column(Column(name='star_formation_rate', data=n.log10(sfr), unit='log10(Msun/yr)'))
t.add_column(Column(name='is_quiescent', data=QU, unit=''))
t.add_column(Column(name='LX_hard', data=n.log10(gal_LX), unit='log10(erg/s)'))
t.add_column(Column(name='mag_abs_r', data=n.zeros_like(mass), unit='mag'))
t.add_column(Column(name='mag_r', data=n.zeros_like(mass), unit='mag'))

t.add_column(Column(name='UM_A_UV'       , data = UM_A_UV   , unit='mag'))   
t.add_column(Column(name='UM_True_SM'    , data = UM_SM     , unit='Msun'))     
t.add_column(Column(name='UM_ICL_mass'   , data = UM_ICL    , unit='Msun'))    
t.add_column(Column(name='UM_True_SFR'   , data = UM_SFR    , unit='Msun/yr'))    
t.add_column(Column(name='UM_Obs_UV'     , data = UM_Obs_UV , unit='mag')) 

# A_UV: UV attenuation (mag)
# SM: True stellar mass (Msun)
# ICL: True intracluster stellar mass (Msun)
# SFR: True star formation rate (Msun/yr)
# Obs_SM: observed stellar mass, including random & systematic errors (Msun)
# Obs_SFR: observed SFR, including random & systematic errors (Msun/yr)
# Obs_UV: Observed UV Magnitude (M_1500 AB)


t.write(path_2_galaxy_file, overwrite=True)
print('done', time.time() - t0, 's')

"""
# STAR FORMATION RATE (valid only for star forming galaxies !)
sfr = n.zeros_like(zz)
# whitaker et al 2012, Eq. 1,2,3.
mean_SFR = (0.70 - 0.13 * zz) * (mass - 10.5) + 0.38 + 1.14 * zz - 0.19 * zz**2.
rds2 = norm.rvs(loc=0, scale=0.34, size=N_halos)
log_sfr = mean_SFR + rds2
print('log sfr', log_sfr[:10], time.time() - t0)
"""

"""
# mass-SFR sequence for the quiescent
sfr_Q = n.zeros_like(zz)


def beta_z(z): return -0.57 * z + 1.43


def alpha_z(z): return 6.32 * z - 16.26


def mean_SFR_Q(mass, z): return mass * beta_z(z) + alpha_z(z)


def scale_z(z): return -0.34 * z + 0.99


rds2 = norm.rvs(loc=0, scale=1., size=len(zz[QU])) * scale_z(zz[QU])

log_sfr_Q = mean_SFR_Q(mass[QU], zz[QU]) + rds2
# change SFR for the quiesent selection
log_sfr[QU] = log_sfr_Q

print('N QU, log SFR[:10]', len(log_sfr[QU]),
      log_sfr[QU][:10], time.time() - t0)
print('N SF, log SFR[:10]', len(log_sfr[SF]),
      log_sfr[SF][:10], time.time() - t0)
"""

"""
# galaxy LF abundance matching to Loveday et al. 2016
ngal, M, phi, Err = n.loadtxt(os.path.join(
    os.environ['GIT_AGN_MOCK'], 'data', 'LF_loveday_2015', 'lf.txt'), unpack=True)
Schechter_M_z = interp1d(M, phi * 0.7**3.)
mrange = n.arange(-24.6, -12.2, 0.01)
total_number = n.cumsum(Schechter_M_z(mrange)) * volume * fraction
print('total_number', total_number)

rds = norm.rvs(loc=0, scale=0.15, size=N_galaxies)
print('N_galaxies', N_galaxies)
vmax_galaxy = 10**(log_vmax + rds)

x_itp = n.hstack((0.,
                  total_number[total_number > 0.1],
                  total_number[total_number > 0.1][-1] * 10))
y_itp = n.hstack((mrange[total_number > 0.1][0],
                  mrange[total_number > 0.1],
                  mrange[total_number > 0.1][-1]))
itp_mags = interp1d(x_itp, y_itp)
mags_2_assign = itp_mags(n.arange(N_galaxies) + 1)
print('mags_2_assign', mags_2_assign, mags_2_assign.shape)
id_sort_vmax = n.argsort(vmax_galaxy)[::-1]  # id=0 the lowest vmax
mag_abs_r = mags_2_assign[id_sort_vmax]
print('mag_abs_r mag_r', mag_abs_r[:10], time.time() - t0)

# compute distance modulus
z_array = n.arange(n.min(zz) - 0.05, n.max(zz) + 0.05, 0.0001)
dist_mods = interp1d(z_array, cosmo.distmod(z_array).value)
mag_r = mag_abs_r + dist_mods(zz)
print('mag_r', mag_r[:10], time.time() - t0)

# mass size relation

# https://arxiv.org/pdf/1411.6355.pdf
# Table 2 and 3. r mag line for sersic selection
# equation 2

# http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
b4 = 7.669
b1 = 1.678
f_14_dev = lambda r12 : gammainc(8, b4*(0.7/r12)**(1./6.))
f_14_exp = lambda r12 : gammainc(2, b1*(0.7/r12)**(1./2.)) # From Raichoor 2017
f_14_test = lambda r12, nn : gammainc(2, b1*(0.7/r12)**(1./nn))

re_dev = lambda M_star : 0.16*(M_star)**(0.1) * (1+M_star/(2.42*10**(10)))**(0.76-0.1)
re_exp = lambda M_star : 0.08*(M_star)**(0.16) * (1+M_star/(17.1*10**(10)))**(0.81-0.16)
# interpolates projected size on the sky
radius_2_arcsec = interp1d(n.arange(0.00001,6.5,0.001), lc_setup.cosmoMD.arcsec_per_kpc_proper( n.arange(0.00001,6.5,0.001) ).value)

radius_dev = re_dev(stellar_mass)* radius_2_arcsec(f1['/sky_position/redshift_R'].value[gal][ok])
radius_exp = re_exp(stellar_mass)* radius_2_arcsec(f1['/sky_position/redshift_R'].value[gal][ok])

radius = radius_exp

dm_dev = -2.5*n.log10(f_14_dev(radius))
dm_exp = -2.5*n.log10(f_14_exp(radius))
# initialize the fibermag with the exponential profile
fiber_mag = rmag + dm_exp
"""

"""
### Option: FIGURES ###
if make_figure:
	import matplotlib
	matplotlib.use('Agg')
	matplotlib.rcParams.update({'font.size': 14})
	import matplotlib.pyplot as p

	fig_dir = os.path.join(
		os.environ['GIT_AGN_MOCK'],
		'figures',
		env,
		'galaxy',
	)
	if os.path.isdir(fig_dir) == False:
		os.system('mkdir -p ' + fig_dir)

	f2 = fits.open(path_2_coordinate_file) # [1].data
	zz = f2[1].data['redshift_R']
	sel = (rds < 1e5 / N_obj)

	# histograms

	fig_out = os.path.join(fig_dir, 'SMHMR_mass_hist_' + baseName + '.png')

	p.figure(1, (6., 5.5))
	p.tight_layout()
	X = mass
	p.hist(X, histtype='step', label='all', rasterized=True, lw=4)
	X = mass[QU]
	p.hist(X, histtype='step', label='QU', rasterized=True, lw=2, ls='dashed')
	X = mass[SF]
	p.hist(X, histtype='step', label='SF', rasterized=True, lw=2)
	p.title(baseName)
	p.xlabel('SMHMR_mass')
	p.ylabel('Counts')
	p.grid()
	p.yscale('log')
	p.legend(frameon=False, loc=0)
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(
		fig_dir,
		'star_formation_rate_hist_' +
		baseName +
		'.png')

	p.figure(1, (6., 5.5))
	p.tight_layout()
	X = log_sfr
	p.hist(X, histtype='step', label='all', rasterized=True, lw=4)
	X = log_sfr[SF]
	p.hist(X, histtype='step', label='SF', rasterized=True, lw=2, ls='dashed')
	X = log_sfr[QU]
	p.hist(X, histtype='step', label='QU', rasterized=True, lw=2)
	p.title(baseName)
	p.xlabel('galaxy_star_formation_rate')
	p.ylabel('Counts')
	p.grid()
	p.yscale('log')
	p.legend(frameon=False, loc=0)
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(fig_dir, 'LX_hard_hist_' + baseName + '.png')

	X = n.log10(gal_LX)

	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.hist(X, histtype='step', rasterized=True, lw=4)
	p.title(baseName)
	p.xlabel('galaxy_LX_hard')
	p.ylabel('Counts')
	p.grid()
	p.yscale('log')
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(fig_dir, 'mag_abs_r_hist_' + baseName + '.png')

	X = mag_abs_r

	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.hist(X, histtype='step', rasterized=True, lw=4)
	p.title(baseName)
	p.xlabel('mag_abs_r')
	p.ylabel('Counts')
	p.grid()
	p.yscale('log')
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(fig_dir, 'mag_r_hist_' + baseName + '.png')

	X = mag_r

	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.hist(X, histtype='step', rasterized=True, lw=4)
	p.title(baseName)
	p.xlabel('galaxy_mag_r')
	p.ylabel('Counts')
	p.grid()
	p.yscale('log')
	p.savefig(fig_out)
	p.clf()

	# 2D plots

	fig_out = os.path.join(fig_dir, 'mag_r_vs_zz_' + baseName + '.png')

	X = zz[sel]
	Y = mag_r[sel]
	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.plot(X, Y, 'k,', rasterized=True)
	p.title(baseName)
	p.xlabel('redshift R')
	p.ylabel('mag_r')
	p.grid()
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(fig_dir, 'mag_abs_r_vs_zz_' + baseName + '.png')

	X = zz[sel]
	Y = mag_abs_r[sel]
	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.plot(X, Y, 'k,', rasterized=True)
	p.title(baseName)
	p.xlabel('redshift R')
	p.ylabel('mag_abs_r')
	p.grid()
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(fig_dir, 'SMHMR_mass_vs_zz_' + baseName + '.png')

	X = zz[sel]
	Y = mass[sel]
	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.plot(X, Y, 'k,', rasterized=True)
	p.title(baseName)
	p.xlabel('redshift R')
	p.ylabel('log10 SMHMR_mass')
	p.grid()
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(
		fig_dir,
		'star_formation_rate_vs_zz_' +
		baseName +
		'.png')

	X = zz[sel]
	Y = log_sfr[sel]
	p.figure(1, (6., 5.5))
	p.tight_layout()
	p.plot(X, Y, 'k+', rasterized=True)
	p.title(baseName)
	p.xlabel('redshift R')
	p.ylabel('log10 star_formation_rate')
	p.grid()
	p.savefig(fig_out)
	p.clf()

	fig_out = os.path.join(
		fig_dir,
		'star_formation_rate_vs_mass_' +
		baseName +
		'.png')

	p.figure(1, (6., 5.5))
	p.tight_layout()
	X = mass[sel]
	Y = log_sfr[sel]
	p.plot(X, Y, 'b+', rasterized=True)
	X = mass[sel & QU]
	Y = log_sfr[sel & QU]
	p.plot(X, Y, 'r+', label='QU', rasterized=True)
	p.title(baseName)
	p.legend(frameon=False, loc=0)
	p.xlabel('log10(mass)')
	p.ylabel('log10 star_formation_rate')
	p.grid()
	p.savefig(fig_out)
	p.clf()
"""