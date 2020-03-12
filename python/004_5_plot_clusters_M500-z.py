"""
What it does
------------

Plots the cluster model

References
----------

Command to run
--------------

python3 004_5_plot_clusters.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then	 work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

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


env =  sys.argv[1]

test_dir = os.path.join(os.environ[env])
path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')
clu = fits.open(path_2_CLU_SAT_catalog)

p_2_training_sample = os.path.join(os.environ['GIT_AGN_MOCK'], 'data/training_CBP')
z_HI, m_HI = n.loadtxt( os.path.join( p_2_training_sample, 'M-z-HIFLUGCS.txt')  , unpack = True)
z_XC, m_XC = n.loadtxt( os.path.join( p_2_training_sample, 'M-z-XCOP.txt')      , unpack = True)
z_XX, m_XX = n.loadtxt( os.path.join( p_2_training_sample, 'M-z-XXL.txt')       , unpack = True)

spt = fits.open( os.path.join( p_2_training_sample, 'bleem_2015.fit'))
fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'clusters' )
if os.path.isdir(fig_dir) == False:
	os.system('mkdir -p ' + fig_dir)

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

nl = lambda sel : len(sel.nonzero()[0])
s1 = clu[1].data['CLUSTER_FX_soft']>9e-14
nl(s1)
# opti0nal red sequence figure
fig_out = os.path.join(fig_dir, 'M500c-zpng')

p.figure(1, (6., 6.))

#p.plot(clu[1].data['redshift_R'], n.log10(clu[1].data['HALO_M500c']), 'k,', label='mock, 13e6', rasterized=True, alpha=0.1)

p.plot(clu[1].data['redshift_R'][s1], n.log10(clu[1].data['HALO_M500c'][s1]), 'm+', label='FX>9e-14, 1e5', rasterized=True, alpha=0.1)

p.plot(spt[1].data['z'], n.log10(spt[1].data['M500c']*1e14),  'kx', label='SPT', rasterized=True)
p.plot(z_HI, m_HI, 'gs', label='HIFLUGCS', rasterized=True)
p.plot(z_XC, m_XC, 'b*', label='X-COP', rasterized=True)
p.plot(z_XX, m_XX, 'ro', label='XXL', rasterized=True)


p.legend(frameon=False)
p.xlabel('redshift')
p.ylabel('M500c')
#p.ylabel('probability distribution function')
p.grid()
# p.ylim((0,1.1))
# p.yscale('log')
p.xscale('log')
p.savefig(fig_out)
p.clf()

nl = lambda sel : len(sel.nonzero()[0])
s1 = clu[1].data['CLUSTER_FX_soft']>2e-14
print(nl(s1))

# opti0nal red sequence figure
fig_out = os.path.join(fig_dir, 'M500c-z-1em14.png')

p.figure(1, (6., 6.))

#p.plot(clu[1].data['redshift_R'], n.log10(clu[1].data['HALO_M500c']), 'k,', label='mock, 13e6', rasterized=True, alpha=0.1)

p.plot(clu[1].data['redshift_R'][s1], n.log10(clu[1].data['HALO_M500c'][s1]), 'm+', label='FX>2e-14, 6e5', rasterized=True, alpha=0.1)

p.plot(spt[1].data['z'], n.log10(spt[1].data['M500c']*1e14),  'kx', label='SPT', rasterized=True)
p.plot(z_HI, m_HI, 'gs', label='HIFLUGCS', rasterized=True)
p.plot(z_XC, m_XC, 'b*', label='X-COP', rasterized=True)
p.plot(z_XX, m_XX, 'ro', label='XXL', rasterized=True)


p.legend(frameon=False)
p.xlabel('redshift')
p.ylabel('M500c')
#p.ylabel('probability distribution function')
p.grid()
# p.ylim((0,1.1))
# p.yscale('log')
p.xscale('log')
p.savefig(fig_out)
p.clf()

