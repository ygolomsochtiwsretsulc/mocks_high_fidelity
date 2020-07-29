#############################################
#############################################
#############################################
#
#
# AGN MODEL
# 
#
#############################################
#############################################
#############################################

# Most recent update of the code march 10

# tabulate K-correction for AGN model
# machine: (pyxspec dependency)
# script: $GIT_AGN_MOCK/python/xspec/agn_tabulate_obscured_fraction.py
# inputs: parameters of the AGN X-ray spectral model
# outputs: $GIT_AGN_MOCK/data/xray_k_correction
#    v3_fraction_observed_A15_Obs_hard_Obs_soft_fscat_002.txt
#    v3_redshift0_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt
#    v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt
#    v3_fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt
# last update: March 4 2020
cd $GIT_AGN_MOCK/python/
python 003_5_agn_tabulate_obscured_fraction.py

# tabulate K-correction foir AGN model
# machine: ds54 (pyxspec dependency)
# script: $GIT_AGN_MOCK/python/003_5_agn_create_xspectra.py
# inputs: parameters of the AGN X-ray spectral model
# outputs: X-ray spectra in the simput format: ${env}/cat_AGN_SIMPUT/agn_Xspectra/NH??.?_Z?.?_N*.fits
# last update: Sept 2019
cd $GIT_AGN_MOCK/python/
python 003_5_agn_create_xspectra.py

###########################################
# full sky AGN model. 
# Full eROSITA simulations.
# script: run_md_agn.py
###########################################

# STEP 1:
# Runs the AGN model
# script: 003_0_agn.py
# inputs: ${env}/fits/all_?.?????.fits, ${env}/fits/all_?.?????_coordinates.fits, ${env}/fits/all_?.?????_galaxy.fits
# outputs: ${env}/fits/all_?.?????_agn.fits

# STEP 2:
# script: 003_1_agn_catalogs.py
# inputs: ${env}/fits/all_?.?????_agn.fits
# outputs:
# 768 files in : ${env}/fits/cat_AGN_?.?????/000???.fit 

# STEP 3: takes about 6 hours
# concatenate all pixels from shells into pixels with all the redshifts
# script: TOPCAT commands in run_md_agn.py
# inputs: 
# N_shells x 768 files in : ${env}/fits/cat_AGN_?.?????/000???.fit 
# outputs: 
# 768 files in : ${env}/cat_AGN_all/000???.fit 

# STEP 4:
# compute XLF, logNlogS, logNlogR
# script: 003_2_agn_compute_XLF_logNlogS_R.py
# inputs: 
# 768 files in : ${env}/cat_AGN_all/000???.fit 
# outputs: 
# three folders
# ${env}/cat_AGN_all/
#   logNlogR
#   logNlogS
#   XLF

# STEP 5:
# plots XLF, logNlogS, logNlogR
# script: 003_3_agn_plot_logNlogS.py
# inputs: 
# three folders
# ${env}/cat_AGN_all/
#   logNlogR
#   logNlogS
#   XLF
# outputs: 
# ${env}/cat_AGN_all/figures

nohup nice -n 19 python run_md_agn.py MD04           > logs/run_agn_MD04_all.log & 
nohup nice -n 19 python run_md_agn.py MD10           > logs/run_agn_MD10_all.log & 
nohup nice -n 19 python run_md_agn.py UNIT_fA1_DIR   > logs/run_agn_UNIT_fA1_DIR_all.log &       
nohup nice -n 19 python run_md_agn.py UNIT_fA2_DIR   > logs/run_agn_UNIT_fA2_DIR_all.log &       
nohup nice -n 19 python run_md_agn.py UNIT_fA1i_DIR  > logs/run_agn_UNIT_fA1i_DIR_all.log & 
nohup nice -n 19 python run_md_agn.py UNIT_fA2i_DIR  > logs/run_agn_UNIT_fA2i_DIR_all.log & 


# create AGN SIMPUT files for eROSITA. 
# Takes 2.5 hours
# Considers a subset of the AGNs in the simulation: FX_soft>2e-17
# script: 003_4_agn_simput.py
# inputs: 768 files in : ${env}/cat_AGN_all/000???.fit 
# outputs: 768 files in : ${env}/cat_AGN_SIMPUT/000???.fit 
# public link (35G)
#  cp -r $MD10/cat_AGN_SIMPUT /data42s/comparat/firefly/mocks/2020-03/MDPL2/
# 
# alternative script, per pixel, more parallel: python 003_4_agn_simput_HPX.py MD10 767
nohup nice -n 19 python 003_4_agn_simput.py MD10  1e-17  > logs/run_simput_MD10_1e17.log & # 
nohup nice -n 19 python 003_4_agn_simput.py MD10  2e-17  > logs/run_simput_MD10_2e17.log & # 
nohup nice -n 19 python 003_4_agn_simput.py MD04  2e-20   > logs/run_simput_MD04.log & # 
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA1_DIR  1e-17 > logs/run_simput_UNIT_fA1_DIR.log &
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA2_DIR  1e-17 > logs/run_simput_UNIT_fA2_DIR.log &
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA1i_DIR 1e-17 > logs/run_simput_UNIT_fA1i_DIR.log &
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA2i_DIR 1e-17 > logs/run_simput_UNIT_fA2i_DIR.log &

#############################################
# Optical/IR properties of AGN (for a subset)
# on he1srv
#############################################

# writes slurm scripts to be executed on he1srv
# script 003_6_agn_magnitudes_WriteScripts.py
# output: cd $GIT_AGN_MOCK/python/003_6_log/*.sh
python 003_6_agn_magnitudes_WriteScripts.py

# compute magnitudes for all AGN simulated with sixte (FX_soft>2e-17)
# Execute the scripts in the directory 003_6_log
# script: 003_6_agn_magnitudes.py
# inputs: 2 x 768 files in : ${env}/cat_AGN_all/000???.fit, ${env}/cat_AGN_SIMPUT/000???.fit 
# outputs: 768 files in : ${env}/cat_AGN-MAG_all/000???.fit 
# e.g.
nohup python3 003_6_agn_magnitudes.py MD10 all 460 > logs/003_6_mag_460.log &  

# compute SDSS-5 / 4MOST AGN / QSO catalogues for BHM and 4MOST consortium surveys S6 and S8
# script: 003_7_4most_catalog.py
# inputs: 768 files in : ${env}/cat_AGN-MAG_all/000???.fit,
# outputs: 768 files in : ${env}/cat_AGN_4MOST/*_000???.fit 
#   global selection :    |g_lat|>10
#      AGN_IR_000???.fit     ( AGN_type=22 | AGN_type=21 ) < 19
#      AGN_DEEP_000???.fit   g_lon>180 & |ecl_lat|>80 & FX > flux_limit_SNR3 - 0.1
#      AGN_WIDE_000???.fit   g_lon>180 & |ecl_lat|<80 & FX > flux_limit_SNR3 - 0.1
#      LyA_000???.fit        AGN_type==11 & AGN_SDSS_r_magnitude < 22.6 & z > 2.2
#      QSO_000???.fit        AGN_type==11 & AGN_SDSS_r_magnitude < 23.2 & 0.9 < z < 2.2
# 
cd $GIT_AGN_MOCK/python/
nohup nice -n 19 python 003_7_4most_catalog.py MD04          > logs/agn_4most_MD04.log &  
nohup nice -n 19 python 003_7_4most_catalog.py MD10          > logs/agn_4most_MD10.log &  
nohup nice -n 19 python 003_7_4most_catalog.py UNIT_fA1_DIR  > logs/agn_4most_UNIT_fA1_DIR.log &  
nohup nice -n 19 python 003_7_4most_catalog.py UNIT_fA2_DIR  > logs/agn_4most_UNIT_fA2_DIR.log &  
nohup nice -n 19 python 003_7_4most_catalog.py UNIT_fA1i_DIR > logs/agn_4most_UNIT_fA1i_DIR.log &  

# Concatenate catalogs 
# inputs: 768 files in : ${env}/cat_AGN_4MOST/*_000???.fit 
# outputs: 
#     cat_AGN-MAG_all ==>> ${env}/AGN_all.fits
#     cat_AGN_4MOST/AGN_DEEP* ==>> ${env}/AGN_DEEP_4MOST.fits
#     cat_AGN_4MOST/AGN_WIDE* ==>> ${env}/AGN_WIDE_4MOST.fits
#     cat_AGN_4MOST/AGN_IR* ==>> ${env}/AGN_IR_4MOST.fits
#     cat_AGN_4MOST/QSO* ==>> ${env}/QSO_4MOST.fits
#     cat_AGN_4MOST/LyA* ==>> ${env}/LyA_4MOST.fits
#   public links
#     /data42s/comparat/firefly/mocks/2020-03/QMOST/S6_4MOST_R-228-234.fit
#     /data42s/comparat/firefly/mocks/2020-03/QMOST/4MOST_S8_AGN.fits
#     /data42s/comparat/firefly/mocks/2020-03/MDPL2/MD10_AGN_all.fits
#     /data42s/comparat/firefly/mocks/2020-03/SMDPL/MD04_AGN_all.fits

cd $GIT_AGN_MOCK/python/

ls $MD10/cat_AGN-MAG_all/*.fit > lists/md10_4most_agn_all.lis
nohup stilts tcat in=@lists/md10_4most_agn_all.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_all.fits          > logs/md10_4most_agn_all_CONCAT.log &

ls $MD10/cat_AGN_4MOST/AGN_DEEP_*.fit > lists/md10_4most_agn_deep.lis
nohup stilts tcat in=@lists/md10_4most_agn_deep.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_DEEP_4MOST.fits  > logs/md10_4most_agn_deep_CONCAT.log &

ls $MD10/cat_AGN_4MOST/AGN_WIDE_*.fit > lists/md10_4most_agn_wide.lis
nohup stilts tcat in=@lists/md10_4most_agn_wide.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_WIDE_4MOST.fits  > logs/md10_4most_agn_wide_CONCAT.log &

ls $MD10/cat_AGN_4MOST/AGN_IR_*.fit > lists/md10_4most_agn_ir.lis
nohup stilts tcat in=@lists/md10_4most_agn_ir.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_IR_4MOST.fits      > logs/md10_4most_agn_ir_CONCAT.log &

ls $MD10/cat_AGN_4MOST/QSO_*.fit > lists/md10_4most_qso.lis
nohup stilts tcat in=@lists/md10_4most_qso.lis ifmt=fits omode=out ofmt=fits out=$MD10/QSO_4MOST.fits            > logs/md10_4most_qso_CONCAT.log &

ls $MD10/cat_AGN_4MOST/LyA_*.fit > lists/md10_4most_lya.lis
nohup stilts tcat in=@lists/md10_4most_lya.lis ifmt=fits omode=out ofmt=fits out=$MD10/LyA_4MOST.fits            > logs/md10_4most_lya_CONCAT.log &

# Verifications on AGN catalogs 
# Computes the number of AGN in each catalogue and present in other catalogues (overlap matrix)
# script: 003_7_shared_agn.py
# inputs: 
#     ${env}/AGN_all.fits
#     ${env}/AGN_DEEP_4MOST.fits
#     ${env}/AGN_WIDE_4MOST.fits
#     ${env}/AGN_IR_4MOST.fits
#     ${env}/QSO_4MOST.fits
#     ${env}/LyA_4MOST.fits
# output: latex tables, https://www.overleaf.com/project/5dd687852c646e0001c1aa58
# nohup nice -n 19 python 003_7_shared_agn.py MD10 > logs/agn_overlap_matrix_MD10.log & 
# nohup nice -n 19 python 003_7_shared_agn.py MD04 > logs/agn_overlap_matrix_MD04.log & 

# Merge the AGN catalogues from MD04 and MD10 into a single 4MOST S6 catalogue
# script: 003_7_merge.py
# inputs: from MD04 and MD10 
#     ${env}/AGN_DEEP_4MOST.fits
#     ${env}/AGN_WIDE_4MOST.fits
#     ${env}/AGN_IR_4MOST.fits
# output: /data42s/comparat/firefly/mocks/2020-03/QMOST/S6_4MOST_R-228-234.fit 
nohup nice -n 19 python 003_7_merge.py > logs/agn_merge_4most_catalogue.log &

# Add footprint bits
# script: 005_4_add_footprint.py
# inputs: catalogue file (any)    
# outputs: updates the files with the footprint bit
nohup python 005_4_add_footprint.py S6_4MOST_R-228-234.fit > logs/4most_footprint_S6_4MOST_R-228-234.log &
nohup python 005_4_add_footprint.py QMOST_S8_QSO.fits > logs/4most_footprint_S8_QSO.log &
nohup python 005_4_add_footprint.py QMOST/QMOST_S8_LyA.fits > logs/4most_footprint_S8_LyA.log &

