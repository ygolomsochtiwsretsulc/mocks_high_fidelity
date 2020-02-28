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
import time
t0 = time.time()

from astropy_healpix import healpy
import sys, os, json

from scipy.interpolate import interp1d
from scipy.stats import norm
from astropy.table import Table, Column
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

cosmoMD = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc, Om0=0.307115)
L_box = 1000.0 / 0.6777

import astropy.io.fits as fits
import numpy as n
print('CREATES FITS FILES with SHAM method')
print('------------------------------------------------')
print('------------------------------------------------')

env = sys.argv[1]  # 'MD04'
HEALPIX_id = int(sys.argv[2])
print(env, HEALPIX_id)
doBG = True
doLRG = True
doFILAMENT = True
doELG = True

#doBG = False
#doLRG = False
#doFILAMENT = False
#doELG = False


area = healpy.nside2pixarea(8, degrees=True)

root_dir = os.path.join(os.environ[env])

dir_2_gal_all = os.path.join(root_dir, "cat_GALAXY_all")

dir_2_agn_all = os.path.join(root_dir, "cat_AGN_all")

dir_2_OUT = os.path.join(root_dir, "cat_SHAM_COSMO")

if os.path.isdir(dir_2_OUT) == False:
	os.system('mkdir -p ' + dir_2_OUT)

# dz = 0.005 is a finer binning than any snapshot steps. The smallet step is 0.011 in dz.
dz = 0.05 

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
BG_density     = 200. # 185138
LRG_density    = 400. # 253764
ELG_density    = 1200.
#QSO_density    = 190.
#QSOLya_density = 50.
#BG_S5_density = 400. # 547756
BG_S5_density_unique = 200. # 268067

# defines the total numbers
NN_s8bg  = get_nz(fn, 's8bg' , BG_density * 1.11  * area * dz, dz)
selection_s8bg = (NN_s8bg>8) & (zmid_all>0.15) & (zmid_all<0.41) 
print('NN_s8bg N/deg2 = ', n.sum(NN_s8bg[selection_s8bg])/area )

NN_s8lrg = get_nz(fn, 's8lrg', LRG_density * 1.14 * area * dz, dz)
selection_s8lrg = (NN_s8lrg>8) & (zmid_all>0.4) & (zmid_all<0.81) 
print('NN_s8lrg N/deg2 = ', n.sum(NN_s8lrg[selection_s8lrg])/area )

NN_s8elg = get_nz(fn, 's8elg', ELG_density * 1.9 * area * dz, dz)
selection_s8elg = (NN_s8elg>8) & (zmid_all>0.6) & (zmid_all<1.11) 
print('NN_s8elg N/deg2 = ', n.sum(NN_s8elg[selection_s8elg])/area )

s5_fn = os.path.join( os.environ['GIT_AGN_MOCK'], 'data', 'cosmo-4most', 'filamentSurvey.nzgmm.json')
s5_mydict  = json_read(s5_fn)
# redshift array
s5_all_z = n.arange(s5_mydict['zmin'], s5_mydict['zmax']+dz, dz)
s5_zmin_all, s5_zmax_all = s5_all_z, s5_all_z+dz
s5_zmid_all = (s5_zmax_all + s5_zmin_all)/2.

NN_s5bg  = get_nz(s5_fn, 'R195' , BG_S5_density_unique * 1.22 * area * dz, dz)
selection_s5bg = (NN_s5bg>8) & (s5_zmid_all>0.05) & (s5_zmid_all<0.51) 
print('NN_s5bg N/deg2 = ', n.sum(NN_s5bg[selection_s5bg])/area )

N_pixels = healpy.nside2npix(8)
def run_mock(HEALPIX_id):
	#HEALPIX_id=359
	# galaxy catalogs
	path_2_gal_all_catalog = os.path.join(dir_2_gal_all, str(HEALPIX_id).zfill(6) + '.fit')
	#path_2_gal_sat_catalog = os.path.join(dir_2_gal_sat, str(HEALPIX_id).zfill(6) + '.fit')

	# agn catalogs
	path_2_agn_all_catalog = os.path.join(dir_2_agn_all, str(HEALPIX_id).zfill(6) + '.fit')
	#path_2_agn_sat_catalog = os.path.join(dir_2_agn_sat, str(HEALPIX_id).zfill(6) + '.fit')

	# output catalog
	path_2_OUT_catalog_BG  = os.path.join(dir_2_OUT, 'BG_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_BG_S5 = os.path.join(dir_2_OUT, 'S5GAL_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_LRG = os.path.join(dir_2_OUT, 'LRG_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_ELG = os.path.join(dir_2_OUT, 'ELG_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_OUT_catalog_QSO = os.path.join(dir_2_OUT, 'QSO_' + str(HEALPIX_id).zfill(6) + '.fit')

	print('=================================================================')
	print(path_2_gal_all_catalog,  path_2_agn_all_catalog)
	print( path_2_OUT_catalog_BG  )
	print( path_2_OUT_catalog_BG_S5  )
	print( path_2_OUT_catalog_LRG )
	print( path_2_OUT_catalog_ELG )
	print( path_2_OUT_catalog_QSO )

	hd_all = fits.open(path_2_gal_all_catalog)[1].data
	N_GAL_all = len(hd_all['ra'])

	####################################
	####################################
	# GALAXIES
	####################################
	####################################

	ra_all  =  hd_all['ra']           
	dec_all =  hd_all['dec']          
	ebv_all =  hd_all['galactic_ebv'] 

	all_mvir = hd_all['HALO_Mvir']                  
	logm     = hd_all['galaxy_stellar_mass']        
	sfr      = hd_all['galaxy_star_formation_rate'] 
	allz     = hd_all['redshift_R']                 
	rds      = norm.rvs(loc=0, scale=0.15, size=len(logm))
	all_vmax = logm + rds
	N_halos = len(all_vmax)
	
	if doBG:
		# LRG1 CASE
		# SHAM with scatter until NZ is filled for cen and sat
		# scatter 0.15
		lrg1_selection = (n.ones(N_halos)==0)
		print("BG S8, time=", time.time()-t0)
		for zmin, zmax, N_lrg1 in zip(zmin_all[(selection_s8bg)], zmax_all[(selection_s8bg)], NN_s8bg[(selection_s8bg)]):
			z_sel = (allz>=zmin)&(allz<zmax)
			all_vmax_sort_id = n.argsort(all_vmax[z_sel])
			min_mass = all_vmax[z_sel][all_vmax_sort_id[-N_lrg1-1]]
			mass_selection = (all_vmax>min_mass)&(z_sel)
			lrg1_selection = (mass_selection) | (lrg1_selection)

		print('N LRG1 S8 selected', len(lrg1_selection.nonzero()[0]),', density=', len(lrg1_selection.nonzero()[0])/area )
		if os.path.isfile(path_2_OUT_catalog_BG):
			os.system("rm " + path_2_OUT_catalog_BG)
		t = Table()
		t['RA'] = Column(ra_all[lrg1_selection], unit='degree', dtype=n.float64)
		t['DEC'] = Column(dec_all[lrg1_selection], unit='degree', dtype=n.float64)
		t['Z'] = Column(allz[lrg1_selection], unit='', dtype=n.float32)
		t['Mstar'] = Column(logm[lrg1_selection], unit='log10(Mass/[Msun])', dtype=n.float32)
		t['SFR'] = Column(sfr[lrg1_selection], unit='log10(SFR/[Msun/yr])', dtype=n.float32)
		t['EBV'] = Column(ebv_all[lrg1_selection], unit='mag', dtype=n.float32)
		t.write(path_2_OUT_catalog_BG)#
		print(path_2_OUT_catalog_BG, 'written', time.time() - t0)

	if doFILAMENT:
		# S5 BG case
		s5_selection = (n.ones(N_halos)==0)
		print("Filament S5", time.time()-t0)
		for zmin, zmax, N_s5lrg1 in zip(s5_zmin_all[selection_s5bg], s5_zmax_all[selection_s5bg], NN_s5bg[selection_s5bg]):
			z_sel = (allz>=zmin)&(allz<zmax)
			all_vmax_sort_id = n.argsort(all_vmax[z_sel])
			min_mass = all_vmax[z_sel][all_vmax_sort_id[-N_s5lrg1-1]]
			mass_selection = (all_vmax>min_mass)&(z_sel)
			s5_selection = (mass_selection) | (s5_selection)

		print('N LRG1 S5 selected', len(s5_selection.nonzero()[0]), ', density=', len(s5_selection.nonzero()[0])/area )
		# write a BG catalogue
		if os.path.isfile(path_2_OUT_catalog_BG_S5):
			os.system("rm " + path_2_OUT_catalog_BG_S5)
		t = Table()
		t['RA'] = Column(ra_all[s5_selection], unit='degree', dtype=n.float64)
		t['DEC'] = Column(dec_all[s5_selection], unit='degree', dtype=n.float64)
		t['Z'] = Column(allz[s5_selection], unit='', dtype=n.float32)
		t['Mstar'] = Column(logm[s5_selection], unit='log10(Mass/[Msun])', dtype=n.float32)
		t['SFR'] = Column(sfr[s5_selection], unit='log10(SFR/[Msun/yr])', dtype=n.float32)
		t['EBV'] = Column(ebv_all[s5_selection], unit='mag', dtype=n.float32)
		t.write(path_2_OUT_catalog_BG_S5)#
		print(path_2_OUT_catalog_BG_S5, 'written', time.time() - t0)

	if doLRG:
		# LRG2 CASE
		lrg2_selection = (n.ones(N_halos)==0)
		print("LRG 2", time.time()-t0)
		for zmin, zmax, N_lrg2 in zip(zmin_all[(selection_s8lrg)], zmax_all[(selection_s8lrg)], NN_s8lrg[(selection_s8lrg)]):
			z_sel = (allz>=zmin)&(allz<zmax)& (lrg1_selection==False) 
			all_vmax_sort_id = n.argsort(all_vmax[z_sel])
			min_mass = all_vmax[z_sel][all_vmax_sort_id[-N_lrg2-1]]
			mass_selection = (all_vmax>min_mass) & (z_sel)
			rds = n.random.rand(len(all_vmax))
			lrg2_selection = (mass_selection) | (lrg2_selection)

		print('N LRG2 selected', len(lrg2_selection.nonzero()[0]), ', density=', len(lrg2_selection.nonzero()[0])/area )
		if os.path.isfile(path_2_OUT_catalog_LRG):
			os.system("rm " + path_2_OUT_catalog_LRG)
		t = Table()
		t['RA'] = Column(ra_all[lrg2_selection], unit='degree', dtype=n.float64)
		t['DEC'] = Column(dec_all[lrg2_selection], unit='degree', dtype=n.float64)
		t['Z'] = Column(allz[lrg2_selection], unit='', dtype=n.float32)
		t['Mstar'] = Column(logm[lrg2_selection], unit='log10(Mass/[Msun])', dtype=n.float32)
		t['SFR'] = Column(sfr[lrg2_selection], unit='log10(SFR/[Msun/yr])', dtype=n.float32)
		t['EBV'] = Column(ebv_all[lrg2_selection], unit='mag', dtype=n.float32)
		t.write(path_2_OUT_catalog_LRG)#
		print(path_2_OUT_catalog_LRG, 'written', time.time() - t0)

	if doELG:
		# ELG parameters
		# ELG select on Mvir
		elg_selection = (n.ones(N_halos)==0)
		mh_mean, mh_scatter = 12.1, 0.3
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

		print('N ELG selected', len(elg_selection.nonzero()[0]),', density=', len(elg_selection.nonzero()[0])/area )
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


#for HEALPIX_id in n.arange(N_pixels):
run_mock(HEALPIX_id)

#for HEALPIX_id in n.arange(N_pixels):
	#print("""nohup nice -n 19  python 005_0_sham_cosmology_catalogs.py MD10 """+str(HEALPIX_id).zfill(3)+""" >  cosmo_4most_run_MD10_"""+str(HEALPIX_id).zfill(3)+""".log &  """   )               
