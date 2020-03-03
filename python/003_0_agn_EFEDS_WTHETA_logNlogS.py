"""
What it does
------------

python 003_0_agn_EFEDS_WTHETA_logNlogS.py

import numpy as n
import os, sys, glob

# about 1e6 full sky randoms :
N_data = 1600 * 100 
area_ratio = 1600 * n.pi/129600.
size = int(N_data * 10 / area_ratio) # 10000 # 40000000

uu = n.random.uniform(size=size)
dec = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
ra = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi

# selection of the region
selection = lambda ra_val, dec_val :  (abs(dec_val)<10) & (ra_val > 130) & (ra_val < 210 )

in_efeds = selection(ra, dec)

dir_2_clustering = os.path.join(os.environ['GIT_2PCF'],'data/randoms')

DATA = n.transpose( [ 
	ra[in_efeds] , 
	dec[in_efeds] ])
n.savetxt(os.path.join(dir_2_clustering, 'randoms.txt'), DATA )

# about 2e6 full sky randoms :
size = int(2500000)

uu = n.random.uniform(size=size)
dec = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
ra = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi

selection = (abs(dec)<10) & (ra > 130) & (ra < 210 )

fraction = len(ra[selection])/len(ra)

full_sky = 129600./n.pi

area =  full_sky * fraction
area = 1585.05
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
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
#import h5py
import numpy as n
print('Creates AGN mock catalogue ')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

env = 'MD10'
path_2_coordinate_file = os.path.join(os.environ[env], 'EFEDS_agn_fsat_10.fits' )

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

z = f1['redshift_R']
lx = f1['LX_hard']
logNH = f1['logNH']
fx = f1['FX_soft']
fx_att = f1['FX_soft_attenuated']
lx_0520 = f1['LX_soft']
logm = f1['SMHMR_mass']
lsar = lx - logm

area = 1585.05

n_agn = len(z)
indexes = n.arange(n_agn)

# TABULATE LOG N LOG S, fx

def get_lognlogs_replicas(fx, area):
	log_f_05_20 = n.log10(fx[fx > 0])
	out = n.histogram(log_f_05_20, bins=n.arange(-20, -8., 0.2))
	# cumulative number density per square degrees
	#x_out = 0.5*(out[1][1:] + out[1][:-1])
	# n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ])
	N_out = n.cumsum(out[0][::-1])[::-1]
	# n.array([n.sum(out[0][ii:]) for ii in range(len(out[0])) ]) / area
	c_out = N_out / area
	c_out_up = (1 + N_out**(-0.5)) * c_out
	c_out_low = (1 - N_out**(-0.5)) * c_out
	c_err = (n.log10(c_out_up) - n.log10(c_out_low)) / 2.
	return N_out  # , c_err

#Xgal = (abs(g_lat)>20)
N_out = get_lognlogs_replicas(fx, area=area)
N_out_att = get_lognlogs_replicas(fx_att, area=area)
N_out_hard = get_lognlogs_replicas(f1['FX_hard'], area=area)
# DATA_X = n.transpose(out)#[(out[1]>0) & (out[2]!=n.inf)]
n.savetxt(path_2_coordinate_file[:-5] + 'logNlogS_soft.ascii'    , N_out)
n.savetxt(path_2_coordinate_file[:-5] + 'logNlogS_soft_att.ascii', N_out_att)
n.savetxt(path_2_coordinate_file[:-5] + 'logNlogS_hard.ascii'    , N_out_hard)

plot_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'eFEDs')

fx_bins = n.arange(-20, -8., 0.2)
x_fx = fx_bins[:-1] + 0.1

def get_lognlogs(ff, area = area):
    outout = n.loadtxt(ff, unpack=True)
    NN = outout / area
    return x_fx, n.log10(NN) # itp


p.figure(1, (6, 6))

# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_AGN_MOCK"],
    'data',
    'logNlogS',
    'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='g', label='G08')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_AGN_MOCK"],
    'data',
    'logNlogS',
    'logNlogS_Merloni_12_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='r', label='M12')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_AGN_MOCK"],
    'data/logNlogS/logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(x_data, n.log10(y_data), lw=3, ls='dotted', color='b', label='M08')

x, y = get_lognlogs(path_2_coordinate_file[:-5] + 'logNlogS_soft.ascii'    )
p.plot(x, y, rasterized=True, lw=1, ls='solid', label='soft')
x, y = get_lognlogs(path_2_coordinate_file[:-5] + 'logNlogS_soft_att.ascii')
p.plot(x, y, rasterized=True, lw=1, ls='solid', label='soft att')
#x, y = get_lognlogs(path_2_coordinate_file[:-5] + 'logNlogS_hard.ascii'    )
#p.plot(x, y, rasterized=True, lw=1, ls='solid', label='hard')

p.xlabel('log(F_X[0.5-2 keV])')
p.ylabel('log(>F_X) [/deg2]')
p.legend(frameon=False, loc=0)
# p.yscale('log')
p.xlim((-17, -11.5))
p.ylim((-2, 4.2))
# p.title('Mocks')
p.grid()
p.savefig(os.path.join(plot_dir, "logN_logS_AGN_soft.png"))
p.clf()


p.figure(1, (6, 6))

# Merloni 2012
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_AGN_MOCK"],
    'data',
    'logNlogS',
    'logNlogS_Merloni_12_AGN_hard.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='r', label='M12')

x, y = get_lognlogs(path_2_coordinate_file[:-5] + 'logNlogS_hard.ascii'    )
p.plot(x, y, rasterized=True, lw=1, ls='solid', label='hard')

p.xlabel('log(F_X[2-10 keV])')
p.ylabel('log(>F_X) [/deg2]')
p.legend(frameon=False, loc=0)
# p.yscale('log')
p.xlim((-17, -11.5))
p.ylim((-2, 4.2))
# p.title('Mocks')
p.grid()
p.savefig(os.path.join(plot_dir, "logN_logS_AGN_hard.png"))
p.clf()
