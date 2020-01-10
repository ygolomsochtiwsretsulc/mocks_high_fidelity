
# concatenate these catalogs 
#cat_AGN-MAG_all
#cat_AGN-MAG_sat
#simput_data = fits.open(path_2_sipmut)
#sat = SRC_ID >= 2e9
#cen = (sat==False)
#indexes_all = (n.arange(N_agn_all) + 1e9 ).astype('int')
#indexes_sat = (n.arange(N_agn_sat) + 2e9 ).astype('int')

cat_AGN-MAG_all/000347.fit 
cat_AGN-MAG_all/000348.fit 
cat_AGN-MAG_all/000349.fit 
cat_AGN-MAG_all/000378.fit 
cat_AGN-MAG_all/000379.fit 
cat_AGN-MAG_all/000442.fit 
cat_AGN-MAG_all/000443.fit 
cat_AGN-MAG_all/000444.fit 
cat_AGN-MAG_all/000445.fit 
cat_AGN-MAG_all/000314.fit 
cat_AGN-MAG_all/000315.fit 
cat_AGN-MAG_all/000316.fit 
cat_AGN-MAG_all/000317.fit 
cat_AGN-MAG_all/000380.fit 
cat_AGN-MAG_all/000381.fit 
cat_AGN-MAG_all/000411.fit 
cat_AGN-MAG_all/000412.fit 
cat_AGN-MAG_all/000413.fit 
cat_AGN-MAG_sat/000347.fit
cat_AGN-MAG_sat/000348.fit
cat_AGN-MAG_sat/000349.fit
cat_AGN-MAG_sat/000378.fit
cat_AGN-MAG_sat/000379.fit
cat_AGN-MAG_sat/000442.fit
cat_AGN-MAG_sat/000443.fit
cat_AGN-MAG_sat/000444.fit
cat_AGN-MAG_sat/000445.fit
cat_AGN-MAG_sat/000314.fit
cat_AGN-MAG_sat/000315.fit
cat_AGN-MAG_sat/000316.fit
cat_AGN-MAG_sat/000317.fit
cat_AGN-MAG_sat/000380.fit
cat_AGN-MAG_sat/000381.fit
cat_AGN-MAG_sat/000411.fit
cat_AGN-MAG_sat/000412.fit
cat_AGN-MAG_sat/000413.fit

$MD10/efeds.list

stilts tcat in=@efeds.list ifmt=fits omode=out ofmt=fits out=$MD10/AGN_eFEDS.fits 


/home/tdwelly/erosita/eFEDS/legacysurvey/dr8/sweep/sweep_eFEDS.fits
/home/comparat/data/legacysurvey/sweep_eFEDS.fits

g = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[1] / dat['decam_mw_transmission'].T[1])
r_mag = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[2] / dat['decam_mw_transmission'].T[2])
z_mag = 22.5 - 2.5 * n.log10(dat['decam_flux'].T[4] / dat['decam_mw_transmission'].T[4])
gr = g - r_mag
rz = r_mag - z_mag
noJunk = (dat['brick_primary']) & (dat['decam_anymask'].T[1]==0) & (dat['decam_anymask'].T[2]==0) & (dat['decam_anymask'].T[4]==0) & (dat['TYCHO2INBLOB']==False)

22.5 - 2.5 *log10( FLUX_G  / MW_TRANSMISSION_G  )
22.5 - 2.5 *log10( FLUX_R  / MW_TRANSMISSION_R  )
22.5 - 2.5 *log10( FLUX_Z  / MW_TRANSMISSION_Z  )
22.5 - 2.5 *log10( FLUX_W1 / MW_TRANSMISSION_W1 )
22.5 - 2.5 *log10( FLUX_W2 / MW_TRANSMISSION_W2 )
         

#from astropy_healpix import healpy
#import sys
#import os
#import time
#from scipy.interpolate import interp1d
#import astropy.io.fits as fits
#import h5py
#import numpy as n
#print('CREATES FITS FILES')
#print('------------------------------------------------')
#print('------------------------------------------------')
#t0 = time.time()
##import astropy.io.fits as fits
## import all pathes

#env = 'MD10' # sys.argv[1]  # 'MD04'
##f_sat = float(sys.argv[2])  # 0.2
##laptop = sys.argv[3]
##print(env, f_sat) #, laptop)

#stilts_cmd = 'stilts'

#root_dir = os.path.join(os.environ[env])

#dir_2_eRO_all = os.path.join(root_dir, "cat_AGN_all")
#dir_2_eRO_sat = os.path.join(root_dir, "cat_AGN_sat")
#dir_2_eRO_mag_all = os.path.join(root_dir, "cat_AGN-MAG_all")
#dir_2_eRO_mag_sat = os.path.join(root_dir, "cat_AGN-MAG_sat")
#dir_2_SMPT    = os.path.join(root_dir, "cat_AGN_SIMPUT")
#dir_2_multiL  = os.path.join(root_dir, "cat_AGN_multiW")

#if os.path.isdir(dir_2_multiL) == False:
    #os.system('mkdir -p ' + dir_2_SMPT)

#HEALPIX_id = 347
##N_pixels = healpy.nside2npix(8)
##for HEALPIX_id in n.arange(N_pixels):
#path_2_eRO_all_catalog = os.path.join(dir_2_eRO_all, str(HEALPIX_id).zfill(6) + '.fit')
#path_2_eRO_sat_catalog = os.path.join(dir_2_eRO_sat, str(HEALPIX_id).zfill(6) + '.fit')

#path_2_eRO_mag_all_catalog = os.path.join(dir_2_eRO_mag_all, str(HEALPIX_id).zfill(6) + '.fit')
#path_2_eRO_mag_sat_catalog = os.path.join(dir_2_eRO_mag_sat, str(HEALPIX_id).zfill(6) + '.fit')

#path_2_SMPT_catalog = os.path.join(
	#dir_2_SMPT, 'SIMPUT_' + str(HEALPIX_id).zfill(6) + '_1024.fit')
#path_2_multiL_catalog = os.path.join(
	#dir_2_multiL, str(HEALPIX_id).zfill(6) + '.fit')
#print('========================         in            ================================')
#print(path_2_eRO_all_catalog, path_2_eRO_sat_catalog, path_2_SMPT_catalog)
#print('========================         out            ================================')
#print(path_2_multiL_catalog)

#hd_all  = fits.open(path_2_eRO_all_catalog)
#hd_sat  = fits.open(path_2_eRO_sat_catalog)
#hd_mag_all  = fits.open(path_2_eRO_mag_all_catalog)
#hd_mag_sat  = fits.open(path_2_eRO_mag_sat_catalog)
#hd_smpt = fits.open(path_2_SMPT_catalog)

#sat = hd_smpt[1].data['SRC_ID'] >= 2e9
#cen = (sat==False)

#id_all = (hd_smpt[1].data['SRC_ID'][cen] - 1e9).astype('int')
#id_sat = (hd_smpt[1].data['SRC_ID'][sat] - 2e9).astype('int')
    
#data_all = hd_all[1].data[id_all]
#data_sat = hd_sat[1].data[id_sat]
