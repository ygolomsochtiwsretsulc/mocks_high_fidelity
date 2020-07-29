cd $GIT_AGN_MOCK/python/

# Computes COORDINATES and basic GALAXY properties
#
# Step 1: computes coordinates
# script: 001_coordinates.py 
# inputs: ${env}/fits/all_?.?????.fits
# outputs: ${env}/fits/all_?.?????_coordinates.fits
#
# Step 2: computes approximate values of stellar mass and star formation rate
# script: 002_0_galaxy.py
# inputs: ${env}/fits/all_?.?????.fits
# outputs: ${env}/fits/all_?.?????_galaxy.fits 
nohup nice -n 19 python run_md_gal.py MD04          > run_gal_MD04_all.log & #  DONE           
nohup nice -n 19 python run_md_gal.py MD10          > logs/run_gal_MD10_all.log & #  DONE           
nohup nice -n 19 python run_md_gal.py MD40          > logs/run_gal_MD40_all.log & #  DONE           
nohup nice -n 19 python run_md_gal.py UNIT_fA1_DIR  > logs/run_gal_UNIT_fA1_DIR_all.log & #  DONE  
nohup nice -n 19 python run_md_gal.py UNIT_fA2_DIR  > logs/run_gal_UNIT_fA2_DIR_all.log & #  DONE 
nohup nice -n 19 python run_md_gal.py UNIT_fA1i_DIR > logs/run_gal_UNIT_fA1i_DIR_all.log & # DONE 
nohup nice -n 19 python run_md_gal.py UNIT_fA2i_DIR > logs/run_gal_UNIT_fA2i_DIR_all.log & # DONE 


# PIXELIZED GALAXY CATALOGS run 
#
# Step 1: pixelizes the catalogue by combining columns from the dark matter halo file, the coordinate and the galaxy files
# script: 002_1_galaxy_catalogs.py 
# inputs: ${env}/fits/all_?.?????.fits, ${env}/fits/all_?.?????_coordinates.fits, ${env}/fits/all_?.?????_galaxy.fits
# outputs:
# 768 files in : ${env}/fits/cat_GALAXY_?.?????/000???.fit 
# 
# Step 2
# Then concatenates the all healpix catalogs into $MD10/cat_GALAXY_all/
# TOPCAT commands
nohup nice -n 19 python run_galaxy_pixelized_catalogs.py MD04          > logs/run_gal_pixelize_MD04_all.log & #             
nohup nice -n 19 python run_galaxy_pixelized_catalogs.py MD10          > logs/run_gal_pixelize_MD10_all.log & #            
nohup nice -n 19 python run_galaxy_pixelized_catalogs.py MD40          > logs/run_gal_pixelize_MD40_all.log & #            
nohup nice -n 19 python run_galaxy_pixelized_catalogs.py UNIT_fA1_DIR  > logs/run_gal_pixelize_UNIT_fA1_DIR_all.log & #      
nohup nice -n 19 python run_galaxy_pixelized_catalogs.py UNIT_fA2_DIR  > logs/run_gal_pixelize_UNIT_fA2_DIR_all.log & #      
nohup nice -n 19 python run_galaxy_pixelized_catalogs.py UNIT_fA1i_DIR > logs/run_gal_pixelize_UNIT_fA1i_DIR_all.log & #     


# Add K-band absolute magnitude to the pixelized galaxy catalogs
# script: 002_2_galaxy_Kmag.py
# inputs: ${env}/cat_GALAXY_all/*.fit
# outputs: updates ${env}/cat_GALAXY_all/*.fit
nohup nice -n 19 python 002_2_galaxy_Kmag.py MD04 > logs/run_gal_Kmag_MD04_all.log & # 
nohup nice -n 19 python 002_2_galaxy_Kmag.py MD10 > logs/run_gal_Kmag_MD10_all.log & # 
nohup nice -n 19 python 002_2_galaxy_Kmag.py MD40 > logs/run_gal_Kmag_MD40_all.log & # 
nohup nice -n 19 python 002_2_galaxy_Kmag.py UNIT_fA1_DIR > logs/run_gal_Kmag_UNIT_fA1_DIR_all.log & # 
nohup nice -n 19 python 002_2_galaxy_Kmag.py UNIT_fA2_DIR > logs/run_gal_Kmag_UNIT_fA2_DIR_all.log & # 
nohup nice -n 19 python 002_2_galaxy_Kmag.py UNIT_fA1i_DIR > logs/run_gal_Kmag_UNIT_fA1i_DIR_all.log & # 
