"""
What it does
------------

Computes and tabulates the logNlogS, logNlogR, 

Creates figures with the tabulate data

References
----------

Command to run
--------------

python3 003_2_agn_compute_XLF_logNlogS_R.py environmentVAR ftyp

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

ftyp: 'all' or 'sat', type of AGN populating a central halo or satellite haloes.

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, matplotlib

"""
import time
t0 = time.time()

import glob, os, sys

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

import astropy.io.fits as fits
from astropy_healpix import healpy
from scipy.interpolate import interp1d

#import h5py
import numpy as n

print('TABULATES XLF and logNlogS')
print('------------------------------------------------')
print('------------------------------------------------')

env = sys.argv[1] # 'MD04'
print(env)

agn_catalog_list = sorted( n.array( glob.glob( os.path.join( os.environ[env], 'cat_AGN_all', '*.fit' ) ) ) )
agn_catalog_list.sort()

fig_dir      = os.path.join( os.environ['GIT_AGN_MOCK'], 'figures', env, 'agn' )
logNlogS_dir = os.path.join( fig_dir, 'logNlogS' )
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)
if os.path.isdir(logNlogS_dir) == False:
    os.system('mkdir -p ' + logNlogS_dir)

area = healpy.nside2pixarea(8, degrees=True)

# logN logS: X-ray flux binning
fx_bins = n.arange(-20, -8., 0.2)
x_fx = fx_bins[:-1] + 0.1

print('tabulate all logNlogS attenuated')
for path_2_eRO_catalog in agn_catalog_list:
	print(path_2_eRO_catalog, time.time()-t0)
	str_healpix_id = os.path.basename(path_2_eRO_catalog)[:-4]
	fx = fits.open(path_2_eRO_catalog)[1].data['FX_soft_attenuated']
	# AGN data
	print(" TABULATE LOG N LOG S, fx ", time.time()-t0)
	def get_lognlogs_replicas(fx, area):
		log_f_05_20 = n.log10(fx[fx > 0])
		out = n.histogram(log_f_05_20, bins = fx_bins )
		N_out = n.cumsum(out[0][::-1])[::-1]
		return N_out  
	N_out = get_lognlogs_replicas(fx, area=area)
	n.savetxt( os.path.join( logNlogS_dir, 'logNlogS_softAtt_' + str_healpix_id + '.ascii'), N_out ) 

print('tabulate all logNlogS')
for path_2_eRO_catalog in agn_catalog_list:
	print(path_2_eRO_catalog, time.time()-t0)
	str_healpix_id = os.path.basename(path_2_eRO_catalog)[:-4]
	fx = fits.open(path_2_eRO_catalog)[1].data['FX_soft']
	# AGN data
	print(" TABULATE LOG N LOG S, fx ", time.time()-t0)
	def get_lognlogs_replicas(fx, area):
		log_f_05_20 = n.log10(fx[fx > 0])
		out = n.histogram(log_f_05_20, bins = fx_bins )
		N_out = n.cumsum(out[0][::-1])[::-1]
		return N_out  
	N_out = get_lognlogs_replicas(fx, area=area)
	n.savetxt( os.path.join( logNlogS_dir, 'logNlogS_soft_' + str_healpix_id + '.ascii'), N_out ) 


def get_lognlogs(hpx_id='000061'):
    ff = os.path.join(logNlogS_dir, 'logNlogS_soft_' + hpx_id + '.ascii')
    outout = n.loadtxt(ff, unpack=True)
    NN = outout / area
    return x_fx, n.log10(NN) 

def get_lognlogs_att(hpx_id='000061'):
    ff = os.path.join(logNlogS_dir, 'logNlogS_softAtt_' + hpx_id + '.ascii')
    outout = n.loadtxt(ff, unpack=True)
    NN = outout / area
    return x_fx, n.log10(NN) 

print('reads all logNlogS and makes the figure')
hpx_ids = n.array([str(el).zfill(6) for el in n.arange(768)])

p.figure(1, (6, 6))
ys=[]
for hpx_id in hpx_ids[:10]:
    x, y = get_lognlogs_att(hpx_id)
    p.plot(x, y, rasterized=True, lw=0.5, ls='solid') # label=hpx_id,
    ys.append(y)

p.plot(x, n.median(n.array(ys),axis=0), lw=2, ls='dashed', label='AGN MD10 Att') 

ys=[]
for hpx_id in hpx_ids[:10]:
    x, y = get_lognlogs(hpx_id)
    p.plot(x, y, rasterized=True, lw=0.5, ls='solid') # label=hpx_id,
    ys.append(y)

p.plot(x, n.median(n.array(ys),axis=0), lw=2, ls='dashed', label='AGN MD10') 
# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
	os.environ["GIT_VS"],
	'data',
	'logNlogS',
	'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='g', label='G08')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
	os.environ["GIT_VS"],
	'data',
	'logNlogS',
	'logNlogS_Merloni_12_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='r', label='M12')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
	os.environ["GIT_VS"],
	'data/logNlogS/logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(x_data, n.log10(y_data), lw=3, ls='dotted', color='b', label='M08')

p.xlabel('log10(F_X[0.5-2 keV])')
p.ylabel('log10(>F_X) [/deg2]')
p.legend(frameon=False, loc=0)
# p.yscale('log')
p.xlim((-19, -11.5))
p.ylim((-2, 5))
# p.title('Mocks')
p.grid()
p.savefig(os.path.join(logNlogS_dir, "logN_logS_AGN.png"))
p.clf()


p.figure(1, (6, 6))

ys=[]
for hpx_id in hpx_i:
    x, y = get_lognlogs(hpx_id)
    #p.plot(x, y, rasterized=True, lw=0.5, ls='solid') # label=hpx_id,
    ys.append(y)

#p.plot(x, n.median(n.array(ys),axis=0), lw=2, ls='dashed', label='AGN MD10') 

ref_line = interp1d(x, n.median(n.array(ys),axis=0))

ys=[]
for hpx_id in hpx_ids:
    x, y = get_lognlogs_att(hpx_id)
    p.plot(x, y-ref_line(x), rasterized=True, lw=0.5, ls='solid') # label=hpx_id,
    ys.append(y)

p.plot(x, n.median(n.array(ys),axis=0)-ref_line(x), lw=2, ls='dashed', label='AGN MD10 Att') 

# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
	os.environ["GIT_VS"],
	'data',
	'logNlogS',
	'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data)-ref_line(n.log10(x_data)), lw=3, ls='dotted', color='g', label='G08')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
	os.environ["GIT_VS"],
	'data',
	'logNlogS',
	'logNlogS_Merloni_12_AGN.data')
x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(n.log10(x_data), n.log10(y_data)-ref_line(n.log10(x_data)), lw=3, ls='dotted', color='r', label='M12')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
	os.environ["GIT_VS"],
	'data/logNlogS/logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = n.loadtxt(path_2_logNlogS_data, unpack=True)
p.plot(x_data, n.log10(y_data)-ref_line(x_data), lw=3, ls='dotted', color='b', label='M08')

p.xlabel('log10(F_X[0.5-2 keV])')
p.ylabel('log10(>F_X)-log10(MD10 AGN)')
p.legend(frameon=False, loc=0)
# p.yscale('log')
p.xlim((-19, -11.5))
p.ylim((-0.4, 0.5))
p.tight_layout()
# p.title('Mocks')
p.grid()
p.savefig(os.path.join(logNlogS_dir, "logN_logS_AGN_ratio.png"))
p.clf()
