"""
What it does
------------

Creates a fits catalog containing cosmology tracers BGx2, LRG, ELG, QSO, LyaQSO for each healpix 768 pixel of 13 deg2 (NSIDE=8)

Command to run
--------------

python3 005_0_sham_cosmology_catalogs.py environmentVAR 

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------
topcat/stilts
import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib

"""
from astropy_healpix import healpy
import sys
import os
import time
from scipy.interpolate import interp1d
from scipy.stats import norm
from astropy.table import Table, Column
import json
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmoMD = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc, Om0=0.307115)
L_box = 1000.0 / 0.6777


import astropy.io.fits as fits
import h5py
import numpy as n
print('CREATES FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits
# import all pathes

env = 'MD10' # sys.argv[1]  # 'MD04'
print(env)

area = healpy.nside2pixarea(8, degrees=True)

stilts_cmd = 'stilts'

root_dir = os.path.join(os.environ[env])

dir_2_gal_all = os.path.join(root_dir, "cat_GALAXY_all")
dir_2_gal_sat = os.path.join(root_dir, "cat_GALAXY_sat")

dir_2_agn_all = os.path.join(root_dir, "cat_AGN_all")
dir_2_agn_sat = os.path.join(root_dir, "cat_AGN_sat")

dir_2_OUT = os.path.join(root_dir, "cat_SHAM_COSMO")

if os.path.isdir(dir_2_OUT) == False:
	os.system('mkdir -p ' + dir_2_OUT)

"""
# decision on the NX depends on the distance between snapshots given in this file for MD10 :
N_snap, Z_snap, A_snap, DC_max, DC_min = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], 'python', 'DM_MDPL2', 'snap_list_with_border.txt' ), unpack=True)
z_array = n.arange(0, 20., 0.001)
dcs = cosmoMD.comoving_distance(z_array)
#z_to_dc = interp1d(z_array, dcs)
dc_to_z = interp1d(dcs, z_array)

Z_min = dc_to_z(DC_min)
Z_max = dc_to_z(DC_max)
id_sort = n.argsort(Z_min)
all_z = Z_min[id_sort]

zmin_all, zmax_all = all_z[:-1], all_z[1:]
zmid_all = (zmax_all + zmin_all)/2.
"""
# dz = 0.005 is a finer binning than any snapshot steps. The smallet step is 0.011 in dz.
# 
dz = 0.005


def json_read(fn):
    with open(fn,'r') as f:
        fixed_json = ''.join(line for line in f if not line.startswith('#'))
        mydict     = json.loads(fixed_json)
    f.close()
    return mydict

def get_nz(fn,sel,ntot,dz):
    # ntot = density * area (total number of targets)
    mydict  = json_read(fn)
    z          = n.arange(mydict['zmin'],mydict['zmax']+dz,dz)
    ps,mus,sds = mydict[sel]['p'],mydict[sel]['mu'],mydict[sel]['sd']
    n1          = n.array([p*norm.pdf(z,mu,sd) for p,mu,sd in zip(ps,mus,sds)])
    n1          = n.sum(n1,axis=0)
    n1         *= ntot
    return n.floor(n1+3).astype('int')

fn = os.path.join( os.environ['GIT_AGN_MOCK'], 'data', 'cosmo-4most', '4most_s8_kidsdr4.nzgmm.json')
mydict  = json_read(fn)
# redshift array
all_z = n.arange(mydict['zmin'],mydict['zmax']+dz, dz)
zmin_all, zmax_all = all_z, all_z+dz
zmid_all = (zmax_all + zmin_all)/2.

# density of tracers per deg-2
BG_density     = 250.
LRG_density    = 400.
ELG_density    = 1200.
QSO_density    = 190.
QSOLya_density = 50.
BG_S5_density = 450. # 547756
BG_S5_density_unique = 200. # 268067

# defines the total numbers
NN_s8bg  = get_nz(fn, 's8bg' , BG_density * 1.27  * area * dz, dz)
selection_s8bg = (NN_s8bg>8) & (zmid_all>0.1) & (zmid_all<0.41) 
print('NN_s8bg', n.sum(NN_s8bg[selection_s8bg])/area )

NN_s8lrg = get_nz(fn, 's8lrg', LRG_density * 1.3 * 1.02 * area * dz, dz)
selection_s8lrg = (NN_s8lrg>8) & (zmid_all>0.4) & (zmid_all<0.81) 
print('NN_s8lrg', n.sum(NN_s8lrg[selection_s8lrg])/area )

NN_s8elg = get_nz(fn, 's8elg', ELG_density * 1.14 * area * dz, dz)
selection_s8elg = (NN_s8elg>8) & (zmid_all>0.6) & (zmid_all<1.11) 
print('NN_s8elg', n.sum(NN_s8elg[selection_s8elg])/area )

NN_s5bg  = get_nz(fn, 's5bg' , BG_S5_density * 1.35 * area * dz, dz)
selection_s5bg = (NN_s5bg>8) & (zmid_all>0.1) & (zmid_all<0.41) 
print('NN_s5bg', n.sum(NN_s5bg[selection_s5bg])/area )

#def get_dn_dv(file_name = os.path.join(os.environ['GIT_EMERGE'], "data/NZ/4FS_scenario2.17Jul2017.CoLoRe.BGhiz.nz")):
  #print( file_name)
  #zmean, dN_dz = n.loadtxt(file_name, unpack=True)
  #dz_all = zmean[1:]-zmean[:-1]
  #dz = dz_all[0] #0.025 # 
  #zmin = zmean-dz/2.
  #zmax = zmean+dz/2.
  #N_per_bin = dN_dz * dz 
  #return zmin, zmax, N_per_bin 


## rehaSH TO SHELL BOUNDARIES
#def rehash(data, z_start=0, z_reach=6):
  #zmin, zmax, N_p_deg2 = data
  #id_min = n.searchsorted(zmax, z_start)
  #dz_low = z_start - zmin[id_min]
  #dz_high = zmax[id_min] - z_start
  #dz_all = zmax[id_min] - zmin[id_min]
  #NP_low = N_p_deg2[id_min]*dz_high/dz_all
  #id_max = n.searchsorted(zmax, z_reach)
  #dz_low = z_reach - zmin[id_max]
  #dz_high = zmax[id_max] - z_reach
  #dz_all = zmax[id_max] - zmin[id_max]
  #NP_high = N_p_deg2[id_max]*dz_low/dz_all
  #bin_sel = (zmax>z_start)&(zmin<z_reach)
  #zmin_o = n.hstack((z_start, zmin[id_min+1:id_max+1]))
  #zmax_o = n.hstack((zmin[id_min+1:id_max+1],z_reach))
  #N_p_deg2_o = n.hstack((NP_low, N_p_deg2[id_min+1:id_max],NP_high))
  #return zmin_o, zmax_o, N_p_deg2_o

#zmin_lrg1, zmax_lrg1, N_lrg1_pdeg2 = rehash(n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], "data/cosmo-4most/bg.nz"), unpack=True))
#N_lrg1_pdeg2 = 2 * N_lrg1_pdeg2 
#zmin_lrg2, zmax_lrg2, N_lrg2_pdeg2 = rehash(n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], "data/cosmo-4most/lrg.nz"), unpack=True) )
#zmin_elg, zmax_elg, N_elg_pdeg2    = rehash(n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], "data/cosmo-4most/elg.nz"), unpack=True)   )
#zmin_qso, zmax_qso, N_qso_pdeg2    = rehash(n.loadtxt(os.path.join(os.environ['GIT_AGN_MOCK'], "data/cosmo-4most/qso.nz"), unpack=True)   )
#print(zmin_qso, zmax_qso, N_qso_pdeg2 )
N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels)[::-1][340:]:
	#HEALPIX_id=359
	# galaxy catalogs
	path_2_gal_all_catalog = os.path.join(dir_2_gal_all, str(HEALPIX_id).zfill(6) + '.fit')
	path_2_gal_sat_catalog = os.path.join(dir_2_gal_sat, str(HEALPIX_id).zfill(6) + '.fit')

	# agn catalogs
	path_2_agn_all_catalog = os.path.join(dir_2_agn_all, str(HEALPIX_id).zfill(6) + '.fit')
	path_2_agn_sat_catalog = os.path.join(dir_2_agn_sat, str(HEALPIX_id).zfill(6) + '.fit')

	# output catalog
	path_2_OUT_catalog_BG  = os.path.join(dir_2_OUT, 'BG_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_BG_S5 = os.path.join(dir_2_OUT, 'S5GAL_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_LRG = os.path.join(dir_2_OUT, 'LRG_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_ELG = os.path.join(dir_2_OUT, 'ELG_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_QSO = os.path.join(dir_2_OUT, 'QSO_' + str(HEALPIX_id).zfill(6) + '.fit')

	print('=================================================================')
	print(path_2_gal_all_catalog, path_2_gal_sat_catalog, path_2_agn_all_catalog, path_2_agn_sat_catalog)
	print( path_2_OUT_catalog_BG  )
	print( path_2_OUT_catalog_BG_S5  )
	print( path_2_OUT_catalog_LRG )
	print( path_2_OUT_catalog_ELG )
	print( path_2_OUT_catalog_QSO )

	hd_all = fits.open(path_2_gal_all_catalog)[1].data
	hd_sat = fits.open(path_2_gal_sat_catalog)[1].data
	N_GAL_all = len(hd_all['ra'])
	N_GAL_sat = len(hd_sat['ra'])

	####################################
	####################################
	####################################
	####################################
	# AGN
	####################################
	####################################
	####################################
	####################################
	"""
	agn_all = fits.open(path_2_agn_all_catalog)[1].data
	agn_sat = fits.open(path_2_agn_sat_catalog)[1].data
	N_AGN_all_i = len(agn_all['ra'])
	N_AGN_sat_i = len(agn_sat['ra'])

	f_sat_agn = 0.1

	N_AGN_sat = int(N_AGN_all_i * f_sat_agn)

	rd_agn_all = n.random.rand(N_AGN_all_i)
	rd_agn_sat = n.random.rand(N_AGN_sat_i)

	sel0 = ((agn_all['AGN_type']==11)|(agn_all['AGN_type']==12))& (agn_all['redshift_S']>0.8) & (agn_all['redshift_S']<3.5) & (agn_all['AGN_SDSS_r_magnitude']<22.7) & ( rd_agn_all < 1 - f_sat_agn )
	sel1 = ((agn_sat['AGN_type']==11)|(agn_sat['AGN_type']==12))& (agn_sat['redshift_S']>0.8) & (agn_sat['redshift_S']<3.5) & (agn_sat['AGN_SDSS_r_magnitude']<22.7) & ( rd_agn_sat < N_AGN_sat*1./N_AGN_sat_i )

	N_agn_c = len(sel0.nonzero()[0])
	N_agn_s = len(sel1.nonzero()[0])

	n_agn = N_agn_c + N_agn_s

	print(N_agn_c, 'AGN t1 centrals')
	print(N_agn_s, 'AGN t1 sat')
	print(n_agn, 'AGN all')
	print( n_agn / area )

	def get_column(key, hdu_all=agn_all,  hdu_sat=agn_sat, sel0=sel0, sel1=sel1):
		return n.hstack(( hdu_all[key][sel0], hdu_sat[key][sel1] ))

	ra_array = get_column('ra')
	dec_array = get_column('dec')
	z_array = get_column('redshift_S')
	mag_r = get_column('AGN_SDSS_r_magnitude')
	ebv = get_column('galactic_ebv')
	#n.histogram(z_array, bins=n.arange(0,6,0.5))            

	# write an AGN catalogue


	if os.path.isfile(path_2_OUT_catalog_QSO):
		os.system("rm " + path_2_OUT_catalog_QSO)

	t = Table()

	t['RA'] = Column(ra_array, unit='degree', dtype=n.float64)
	t['DEC'] = Column(dec_array, unit='degree', dtype=n.float64)
	t['Z'] = Column(z_array, unit='', dtype=n.float32)
	t['MAG'] = Column(mag_r, unit='mag', dtype=n.float32)
	t['EBV'] = Column(ebv, unit='mag', dtype=n.float32)

	t.write(path_2_OUT_catalog_QSO)#
	print(path_2_OUT_catalog_QSO, 'written', time.time() - t0)
	"""

	####################################
	####################################
	####################################
	####################################
	# GALAXIES
	####################################
	####################################
	####################################
	####################################

	ra_all = n.hstack(( hd_all['ra'], hd_sat['ra'] ))
	dec_all = n.hstack(( hd_all['dec'], hd_sat['dec'] ))
	ebv_all = n.hstack(( hd_all['galactic_ebv'], hd_sat['galactic_ebv'] ))

	all_mvir = n.hstack(( hd_all['HALO_Mvir'], hd_sat['HALO_Mvir'] ))
	logm = n.hstack(( hd_all['galaxy_stellar_mass'], hd_sat['galaxy_stellar_mass'] ))
	sfr = n.hstack(( hd_all['galaxy_star_formation_rate'], hd_sat['galaxy_star_formation_rate'] ))
	allz = n.hstack(( hd_all['redshift_R'], hd_sat['redshift_R'] ))
	rds = norm.rvs(loc=0, scale=0.15, size=len(logm))
	all_vmax = logm + rds
	N_halos = len(all_vmax)

	# ELG parameters
	# ELG select on Mvir
	elg_selection = (n.ones(N_halos)==0)
	mh_mean, mh_scatter = 12.1, 0.3
	print("==============")
	print("==============")
	print("ELG", time.time()-t0)
	for zmin, zmax, N_elg in zip(zmin_all[selection_s8elg], zmax_all[selection_s8elg], NN_s8elg[selection_s8elg]):
		z_sel = (allz>=zmin)&(allz<zmax)
		mh_bins = n.arange(mh_mean -2*mh_scatter, mh_mean +2*mh_scatter+0.05, 0.05)
		mh_bins_pos = 0.5*(mh_bins[1:]+mh_bins[:-1])
		proba = lambda x : norm.pdf(x, loc=mh_mean,scale=mh_scatter)
		proba_norm = proba(mh_bins_pos).sum()
		N_2_select_per_bin = (N_elg*proba(mh_bins_pos)/proba_norm).astype('int')
		for id_bin in range(len(mh_bins)-1):
			id_in_bin =(z_sel) & (all_mvir > mh_bins[id_bin]) &( all_mvir < mh_bins[id_bin+1]) 
			N_avail = len(id_in_bin.nonzero()[0])
			rds = n.random.rand(len(all_mvir))
			bin_selection = (id_in_bin)&(rds < N_2_select_per_bin[id_bin]*1./N_avail)
			elg_selection = (bin_selection)|(elg_selection)

	print('ELG selected', len(elg_selection.nonzero()[0]), len(elg_selection.nonzero()[0])/area )
	# write an ELG catalogue
	if os.path.isfile(path_2_OUT_catalog_ELG):
		os.system("rm " + path_2_OUT_catalog_ELG)
	t = Table()
	t['RA'] = Column(ra_all[elg_selection], unit='degree', dtype=n.float64)
	t['DEC'] = Column(dec_all[elg_selection], unit='degree', dtype=n.float64)
	t['Z'] = Column(allz[elg_selection], unit='', dtype=n.float32)
	t['Mstar'] = Column(logm[elg_selection], unit='log10(Mass/[Msun])', dtype=n.float32)
	t['SFR'] = Column(sfr[elg_selection], unit='log10(SFR/[Msun/yr])', dtype=n.float32)
	t['EBV'] = Column(ebv_all[elg_selection], unit='mag', dtype=n.float32)
	t.write(path_2_OUT_catalog_ELG)#
	print(path_2_OUT_catalog_ELG, 'written', time.time() - t0)

