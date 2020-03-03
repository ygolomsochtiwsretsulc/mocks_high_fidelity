"""
What it does
------------

Plots the cluster model, richness vs mass


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
from scipy.optimize import curve_fit
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import numpy as n
print('Create file with galaxies around clusters')
print('=> Abundance matching for magnitudes')
print('=> Red sequence colors')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


env = sys.argv[1]
#delta_crit = '200c' # sys.argv[2]
#print(env, delta_crit)
test_dir = os.path.join(os.environ[env])
fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'clusters', 'galaxies')
all_catalogs = n.array(glob.glob(os.path.join(fig_dir, 'richness_mass_??????.txt')))
all_aexp = n.array([ os.path.basename(el[:-4]).split('_')[-1] for el in all_catalogs ])
all_aexp.sort()

all_aexp_float = n.array([ float(el[0]+'.'+el[1:]) for el in all_aexp ])
all_redshift = n.array([ 1./el - 1. for el in all_aexp_float ])

def mass_2_richness(M200c, redshift): 
	"""
	M200c to richness conversion using the scaling relation 
	Table 2, second line of Capasso et al. 2019 (1812.0609
	"""
	return 39.8 * (M200c/3e14)**(0.98) * ((1+redshift)/(1+0.18))**(-1.08)
data = []
for path_2_cat, z_i in zip(all_catalogs, all_redshift):
	p_i = n.loadtxt(path_2_cat)
	data.append( n.hstack(( z_i, p_i)) )

DATA = n.transpose(data)
data_z = DATA[0]
data_a0 = DATA[1]
data_a1 = DATA[2]


fig_out = os.path.join(fig_dir, 'richness_mass_evolution_1e14.png')
coeff_out = os.path.join(fig_dir, 'richness_mass_evolution_1e14.txt')
n.savetxt(coeff_out, DATA)


p.figure(1, (6., 5.5))
p.axes([0.15, 0.15, 0.8, 0.8])
p.plot(n.log10(1+data_z), data_a0, 'r+', label='a0')
#p.plot(n.log10(1+data_z), data_a1+10, 'b+', label='a1')
#p.title(r'$\bar{z}$=' + str(n.round(z_snap, 3)))
p.xlabel(r'$\log_{10}(1+z)$')
p.ylabel('a0')#, a1+10')
#p.yscale('log')
#p.xscale('log')
p.legend(loc=0, frameon=False)
p.grid()
p.savefig(fig_out)
p.clf()


