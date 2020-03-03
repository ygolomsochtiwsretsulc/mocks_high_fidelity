"""
What it does
------------

Take the complete cluster galaxy satellite FILE:  $ENV_eRO_CLU_SAT.fit
Re-assigns the is_quiescent flag as a function of radius as well as star formation rates.

It increases the number of quiescent galaxies in the vicinity of clusters

References
----------

 * Hennig et al. 2017, https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.4015H
 * Biffi et al. Dolag et al. Magneticum, https://ui.adsabs.harvard.edu/#abs/2018MNRAS.tmp.2317B

Command to run
--------------

python3 004_3_cluster_red_galaxies.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, numpy, astropy


"""

from scipy.optimize import newton
import sys, os, time
from scipy.stats import norm
from scipy.interpolate import interp1d 
#import astropy.io.fits as fits
import numpy as n
from scipy.special import erf
from astropy.table import Table, Column
print('Adjusts red sequence of galaxies around clusters')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

env = sys.argv[1]
baseName = sys.argv[2]
delta_crit = sys.argv[3]
#env="MD10" 
#baseName="all_0.89510"
#delta_crit = '200c'
#delta_crit = '500c'
#delta_crit = 'vir'
#delta_crit = '2rvir'

z_snap = 1./float(baseName.split('_')[1])-1.
aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
print(env, baseName, delta_crit)
test_dir = os.path.join(os.environ[env])
path_2_CLU_SAT_catalog = os.path.join(test_dir, 'fits', baseName + '_galaxiesAroundClusters.fit')


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


# INCREASE QUIESCENT FRACTIONS !

hdu_clu = Table.read(path_2_CLU_SAT_catalog)
x = ((hdu_clu['x']-hdu_clu['HOST_HALO_x'])**2. + (hdu_clu['y']-hdu_clu['HOST_HALO_y'])**2. + (hdu_clu['z']-hdu_clu['HOST_HALO_z'])**2.)**0.5
is_quiescent = hdu_clu['is_quiescent']
zr_CLU = hdu_clu['redshift_R']
mass = hdu_clu['SMHMR_mass']
log_sfr = hdu_clu['star_formation_rate']


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



def frac_old(x, z_cluster): return (erf((-n.log10(x) + 0.1)/0.6)+0.9)*0.38 * (1+z_cluster)**(-0.65)+0.22 

f_red = frac_old(x, zr_CLU)

rds = n.random.random(len(x))

is_red = (rds < f_red)
is_1rvir = (x <= RADIUS/1000.)

print(len(is_red[is_1rvir].nonzero()[0]))
print(len(is_quiescent[is_1rvir].nonzero()[0]))

N_red_1rvir_needed = len(is_red[is_1rvir].nonzero()[0])
N_quiescent_1rvir_already_there = len(is_quiescent[is_1rvir].nonzero()[0])

N_red_1rvir_to_assign = N_red_1rvir_needed - N_quiescent_1rvir_already_there

# find them among the non_quiescent
SF = (is_quiescent == False) & (is_1rvir)
rds2 = n.random.random(len(x[SF]))
is_1rvir2 = (x[SF] < RADIUS[SF]/1000.)
f_red2 = frac_old(x[SF], zr_CLU[SF])


def func(DELTA_X):
    change_2_red = (rds2 < f_red2 - DELTA_X)
    N_new_red = len(change_2_red[is_1rvir2].nonzero()[0])
    # print(N_new_red)
    return N_red_1rvir_to_assign - N_new_red  # , change_2_red


VALS = interp1d(n.array([ func(xxx) for xxx in n.arange(0,1,0.001) ]), n.arange(0,1,0.001) )
#print(VALS(0))
change_2_red = (rds2 < f_red2 - VALS(0))[(x[SF] < RADIUS[SF]/1000.)]

is_quiescent[(SF) & (x < RADIUS/1000.)] = change_2_red

# change the SFR values
# mass-SFR sequence for the quiescent
sfr_Q = n.zeros_like(zr_CLU)


def beta_z(z): return -0.57 * z + 1.43


def alpha_z(z): return 6.32 * z - 16.26


def mean_SFR_Q(mass, z): return mass * beta_z(z) + alpha_z(z)


def scale_z(z): return -0.34 * z + 0.99


rds3 = norm.rvs(loc=0, scale=1., size=len(
    zr_CLU[(SF) & (x < 1)])) * scale_z(zr_CLU[(SF) & (x < 1)])
log_sfr_Q = mean_SFR_Q(mass[(SF) & (x < 1)], zr_CLU[(SF) & (x < 1)]) + rds3
# change SFR for the quiesent selection
log_sfr[(SF) & (x < 1)] = log_sfr_Q


# update the file SFR and is_quiescent columns
hdu_clu['is_quiescent'] = is_quiescent
hdu_clu['galaxy_star_formation_rate'] = log_sfr


hdu_clu.write(path_2_CLU_SAT_catalog, overwrite=True) 
