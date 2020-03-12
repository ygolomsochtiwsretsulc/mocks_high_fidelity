"""
What it does
------------

Computes the AGN model from Comparat et al. 2019.

Re-samples to a given f_sat

python 003_0_agn_EFEDS_WTHETA.py 2
python 003_0_agn_EFEDS_WTHETA.py 4
python 003_0_agn_EFEDS_WTHETA.py 6
python 003_0_agn_EFEDS_WTHETA.py 8
python 003_0_agn_EFEDS_WTHETA.py 10
python 003_0_agn_EFEDS_WTHETA.py 12
python 003_0_agn_EFEDS_WTHETA.py 14
python 003_0_agn_EFEDS_WTHETA.py 16
python 003_0_agn_EFEDS_WTHETA.py 20

"""
import sys
import os
import time
import extinction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from astropy.table import Table, Column
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
#import h5py
import numpy as n
print('Creates AGN mock catalogue ')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


env = sys.argv[1] # 'MD10'
f_sat_pc = int(sys.argv[2]) # 10 # percent
f_sat = f_sat_pc / 100.
print('f sat', f_sat)
test_dir = os.path.join(os.environ[env])

path_2_coordinate_file = os.path.join(test_dir, 'EFEDS_all.fits')
path_2_agn_file = os.path.join(test_dir, 'EFEDS_agn_fsat_'+str(f_sat_pc)+'.fits')

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

if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR" or env == "MD40" or env == "MD10":
    scatter_0 = 1.4
if env == "MD04":
    scatter_0 = 1.0

print('opens coordinate file ', time.time() - t0)
f1 = Table.read( path_2_coordinate_file )
zz_1 = f1['redshift_R']
cen = (f1['pid']==-1)
sat = (cen==False)
#high_z = (zz_1>1.5)
N_galaxies = len(zz_1)
N_galaxies_cen = len(zz_1[cen])
N_galaxies_sat = len(zz_1[sat])

f_duty = interp1d( n.array([0., 0.75, 1.5, 3.5, 10.1]), n.array([0.1, 0.2, 0.3, 0.3, 0.3]))
f_duty_realization = f_duty(zz_1)
active = (n.random.random(size=N_galaxies) <= f_duty_realization)
#active[high_z][:]=True
N_agn = len(zz_1[active])

print('native N, cen, sat, f_sat', N_galaxies, N_galaxies_cen, N_galaxies_sat, N_galaxies_sat*1./N_galaxies_cen)
print('N AGN', N_agn)

native_f_sat = N_galaxies_sat*1./N_galaxies_cen 

N_galaxies_sat/N_galaxies_cen

rds = n.random.random(N_galaxies)
if f_sat > native_f_sat:
	# downsample centrals
	print('downsamples centrals')
	N_cen_goal = N_galaxies_sat / f_sat	
	sel_cen = (rds < N_cen_goal / N_galaxies_cen)
	all_cen = (cen)&(sel_cen)
	keep = (sat)|(all_cen)
	
if f_sat <= native_f_sat:
	# downsample sat
	print('downsamples sat')
	N_sat_goal = N_galaxies_cen * f_sat	
	sel_sat = (rds<N_sat_goal/N_galaxies_sat)
	all_sat = (sat)&(sel_sat)
	keep = (cen)|(all_sat)

f2 = f1[keep]

zz = f2['redshift_R']
dL_cm = f2['dL']
galactic_NH = f2['nH']
galactic_ebv = f2['ebv']
mass = f2['SMHMR_mass']  # log of the stellar mass
#high_z = (zz>1.5)

cen = (f2['pid']==-1)
sat = (cen==False)

N_galaxies = len(zz)
N_galaxies_cen = len(zz[cen])
N_galaxies_sat = len(zz[sat])


f_duty = interp1d( n.array([0., 0.75, 1.5, 3.5, 10.1]), n.array([0.1, 0.2, 0.3, 0.3, 0.3]))
f_duty_realization = f_duty(zz)
active = (n.random.random(size=N_galaxies) <= f_duty_realization)
#active[high_z][:]=True
N_agn = len(zz[active])

print('obtained N, cen, sat, f_sat', N_galaxies, N_galaxies_cen, N_galaxies_sat, N_galaxies_sat*1./N_galaxies_cen)
print('N AGN', N_agn)

f3 = f2[active]

tt = Table(f3)
tt.write(path_2_agn_file, overwrite=True)

