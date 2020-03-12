"""
What it does
------------

Add K band magnitude to the GALAXY catalogues

Command to run
--------------

python3 002_2_galaxy_Kmag.py environmentVAR 

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/fits/

Dependencies
------------
import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib

"""
import time
t0=time.time()

#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 14})
#import matplotlib.pyplot as p

from astropy_healpix import healpy
import sys, os

from scipy.interpolate import interp1d
from scipy.stats import norm
from astropy.table import Table, Column
from scipy.optimize import curve_fit
import linmix
import astropy.io.fits as fits
#import h5py
import numpy as n
print('CREATES FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits
# import all pathes

env = sys.argv[1]  # 'MD04'
print(env)

root_dir = os.path.join(os.environ[env])
plotDir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'logMs-MK-fits')
catalogue_dir = os.path.join(os.environ[env], 'cat_GALAXY_all')

"""
# first, fits the K-band stellr mass relation as a function of redshift
path_2_COSMOS = os.path.join(os.environ['HOME'],'data/COSMOS/photoz_vers2.0_010312.fits')

t = Table.read(path_2_COSMOS)
good = (t['photoz']>0. )&( t['photoz']< 2. ) & ( t['MK']<0 )&( t['MK']>-40 )&( t['mass_med']<14 )&( t['mass_med']>8 )


dz = 0.1
zall = n.arange(0., 2.+dz, dz)
params, covar = [], []
params_lm = []
resid = []
for zmin, zmax in zip(zall[:-1], zall[1:]):
	sel = (good) & (t['photoz']>zmin) & (t['photoz']<=zmax) &( t['mass_med']>8+zmin/2 )
	#fun = lambda x, a, b : a*x+b
	# fixing the slope to avoid jumps with redshift correlated to the second parameters
	fun = lambda x, b : -2.15 * x + b
	x_data, y_data = t['mass_med'][sel], t['MK'][sel]
	popt, pcov= curve_fit(fun, x_data, y_data, p0=(0.))#, sigma=y_err)
	params.append(popt)
	covar.append(pcov)
	#n.random.seed(2)
	#lm = linmix.LinMix(x_data, y_data, K=2)
	#lm.run_mcmc(miniter=200, maxiter=600)
	#params_lm.append([n.median(lm.chain['beta']), n.median(lm.chain['alpha'])])
	#print(n.round(zmin,1),n.round(zmax,1),'curve_fit',n.round(popt,2), n.round(pcov[0][0],4), n.round(pcov[1][1],4))#,'linmix',params_lm[-1])# = 0.05
	print(n.round(zmin,1),n.round(zmax,1),'curve_fit',n.round(popt,2), n.round(pcov[0][0],4))
	# figure
	xs = n.arange(n.min(x_data), n.max(x_data), (n.max(x_data)-n.min(x_data))/20.)
	p.figure(figsize=(6,6))
	p.plot(x_data, y_data, 'g+', alpha=0.5, label='COSMOS, N='+str(len(x_data)))
	#p.errorbar(x_data, y_data, xerr=x_err, yerr=y_err, ls=' ', alpha=0.5)
	#for i in range(0, len(lm.chain), 25):
		#
		#ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
		#p.plot(xs, ys, color='r', alpha=0.02)
	#ys = n.median(lm.chain['alpha']) + xs * n.median(lm.chain['beta'])
	#p.plot(xs, ys, color='k', label=r'$MK=\log_{10}(M_s)\times$'+str(n.round(n.median(lm.chain['beta']),2)) +'+'+str(n.round(n.median(lm.chain['alpha']),2)))
	#p.plot(xs, fun(xs, popt[0], popt[1]), 'b--', label=r'$MK=\log_{10}(M_s)\times$'+str(n.round(popt[0],2)) +'+'+str(n.round(popt[1],2)))
	p.plot(xs, fun(xs, popt[0]), 'b--', label=r'$MK=-2.15\log_{10}(M_s)\times$'+'+'+str(n.round(popt[1],2)))
	p.xlabel(r'$\log_{10}(M_s/M_\odot)$')
	p.ylabel(r'MK')
	p.legend(frameon=False, loc=3)
	p.ylim((-26,-16))
	p.xlim((8., 12.))
	p.grid()
	p.title(str(n.round(zmin,2))+r'$<z<$'+str(n.round(zmax,2)) )#+ r", $F_X>1\times10^{-17}$ ")
	p.savefig(os.path.join( plotDir, "bayesian_mass_Kmag_"+str(n.round(zmin,2))+".png") )
	p.clf()
	# TREATMENT OF THE RESIDUALS
	# data
	dy_B = y_data - fun(x_data, popt[0])
	#y_model_B = n.median(lm.chain['alpha']) + x_data * n.median(lm.chain['beta'])
	#dy_B = y_data-y_model_B
	#fit
	fun = lambda x, loc, scale : norm.cdf(x, loc=loc, scale=scale) 
	# figure
	bins=n.arange(-3,3,0.1)
	xbins=n.arange(-3,3,0.1)[1:]-0.05
	p.figure(2, (6,6))
	n_t2 = p.hist(dy_B , bins=bins, normed=True, cumulative=True, histtype='step', label='residuals' )[0]
	pt2, ct2 = curve_fit(fun, xbins, n_t2, p0=(0,1))
	p.plot(xbins, norm.cdf(xbins, loc=pt2[0], scale=pt2[1]), label='Gaussian N('+str(n.round(pt2[0],2))+', '+str(n.round(pt2[1],2))+')', ls='dashed')
	p.ylabel('counts')
	p.xlabel(r'delta mag')
	p.grid()
	p.legend(frameon=False, loc=0)
	p.title(str(n.round(zmin,2))+r'$<z<$'+str(n.round(zmax,2)))
	p.savefig(os.path.join( plotDir, "hist_residual_mass_Kmag_"+str(n.round(zmin,2))+".png"))
	p.clf()
	print('residuals', pt2)
	resid.append( [pt2,ct2] )
	#return popt, pcov, n.median(lm.chain['beta']), n.median(lm.chain['alpha']), n.std(lm.chain['beta']), n.std(lm.chain['alpha']), residuals

scatter = 0.5
#p_a, p_b = n.transpose(params)
#OUT = n.transpose([zall[:-1], zall[1:],  p_a, p_b])
p_a = n.transpose(params)
OUT = n.transpose([zall[:-1], zall[1:],  p_a[0]])
n.savetxt(os.path.join( plotDir, 'look-up-table.txt'), OUT)
#sys.exit()
"""

# Now assigns K band absolute magnitudes :
#z0,z1,a0,b0 = n.loadtxt(os.path.join( plotDir, 'look-up-table.txt'), unpack=True)
#fun = lambda x, a, b : a*x+b
# 1 param to have a smoother evolution with redshift
z0,z1,b0 = n.loadtxt(os.path.join( plotDir, 'look-up-table.txt'), unpack=True)
fun = lambda x, b : -2.15*x+b

def add_kmag(HEALPIX_id):
	#HEALPIX_id=359
	# catalogs
	path_2_GAL    = os.path.join(catalogue_dir, str(HEALPIX_id).zfill(6) + '.fit')

	t_gal = Table.read(path_2_GAL)
	if 'K_mag_abs' not in t_gal.columns:
		t_gal.add_column(Column(name='K_mag_abs', data=n.zeros_like(t_gal['galaxy_stellar_mass'].data), unit='mag'))

	for zmin, zmax, p_b in zip(z0, z1, b0):
		#print(zmin, zmax, p_b)
		s_gal  = (t_gal['redshift_R']>=zmin)  & (t_gal['redshift_R']<=zmax) 
		mag = lambda x : fun( x, p_b)
		t_gal ['K_mag_abs'][s_gal]  = mag(t_gal ['galaxy_stellar_mass'][s_gal])

	t_gal ['K_mag_abs']+=norm.rvs(loc=0, scale=0.15, size=len(t_gal['redshift_R']))
	t_gal.write  (path_2_GAL   , overwrite=True)

N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels):
	print(HEALPIX_id, time.time()-t0)
	add_kmag(HEALPIX_id)
