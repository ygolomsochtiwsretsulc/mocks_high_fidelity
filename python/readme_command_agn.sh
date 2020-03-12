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

# tabulate K-correction for AGN model
# machine: ds54 (pyxspec dependency)
# script: $GIT_AGN_MOCK/python/xspec/agn_tabulate_obscured_fraction.py
# inputs: parameters of the AGN X-ray spectral model
# outputs: $GIT_AGN_MOCK/data/xray_k_correction
#    v3_fraction_observed_A15_Obs_hard_Obs_soft_fscat_002.txt
#    v3_redshift0_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt
#    v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt
#    v3_fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt
cd $GIT_AGN_MOCK/python/xspec/
python agn_tabulate_obscured_fraction.py

# tabulate K-correction foir AGN model
# machine: ds54 (pyxspec dependency)
# script: $GIT_AGN_MOCK/python/003_5_agn_create_xspectra.py
# inputs: parameters of the AGN X-ray spectral model
# outputs: X-ray spectra in the simput format: ${env}/cat_AGN_SIMPUT/agn_Xspectra/NH??.?_Z?.?_N*.fits
cd $GIT_AGN_MOCK/python/
python 003_5_agn_create_xspectra.py



###########################################
# small sky area AGN model
# AGN clustering in eFEDS. 
###########################################
# Light cone for redshift < 6 & |dec|<10 & 130 < ra < 210 (Small area around eFEDs)
# It keeps only columns to compute the AGN model.
# script: 007_0_EFEDS_WTHETA.py
# inputs: ${env}/fits/all_?.?????.fits, ${env}/fits/all_?.?????_coordinates.fits, ${env}/fits/all_?.?????_galaxy.fits
# outputs: ${env}/fits/all_?.?????_EFEDSwtheta.fits
nohup nice -n 19  python run_eFEDs_AGN.py MD04            >  EFEDS_AGN_MD04.log &                    
nohup nice -n 19  python run_eFEDs_AGN.py MD10            >  EFEDS_AGN_MD10.log &                    
nohup nice -n 19  python run_eFEDs_AGN.py MD40            >  EFEDS_AGN_MD40.log &                    
nohup nice -n 19  python run_eFEDs_AGN.py UNIT_fA1_DIR    >  EFEDS_AGN_UNIT_fA1_DIR.log &                    
nohup nice -n 19  python run_eFEDs_AGN.py UNIT_fA2_DIR    >  EFEDS_AGN_UNIT_fA2_DIR.log &                    
nohup nice -n 19  python run_eFEDs_AGN.py UNIT_fA1i_DIR   >  EFEDS_AGN_UNIT_fA1i_DIR.log &                    

# TOPCAT command to concatenate the eFEDS light cone into a single file
# then the AGN model will be applied and clustering predicted 
# inputs: ${env}/fits/all_?.?????_EFEDSwtheta.fits
# output: ${env}/EFEDS_all.fits
ls $MD10/fits/*_EFEDSwtheta.fits > md10_efeds.lis
stilts tcat in=@md10_efeds.lis ifmt=fits omode=out ofmt=fits out=$MD10/EFEDS_all.fits
rm md10_efeds.lis

# creates eFEDs mocks with different satellites fractions
# Applies the duty cycle
# script: 003_0_agn_EFEDS_WTHETA.py
# input: ${env}/EFEDS_all.fits
# output: ${env}/EFEDS_agn_fsat_*.fits
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 2   >  EFEDS_AGN_MD10_fsat2_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 4   >  EFEDS_AGN_MD10_fsat4_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 6   >  EFEDS_AGN_MD10_fsat6_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 8   >  EFEDS_AGN_MD10_fsat8_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 12  >  EFEDS_AGN_MD10_fsat12_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 14  >  EFEDS_AGN_MD10_fsat14_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 16  >  EFEDS_AGN_MD10_fsat16_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 20  >  EFEDS_AGN_MD10_fsat20_DC.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA.py MD10 10  >  EFEDS_AGN_MD10_fsat10_DC.log &

# Runs the AGN model
# script: 003_0_agn_EFEDS_WTHETA_model.py
# input: ${env}/EFEDS_agn_fsat_*.fits
# output: ${env}/EFEDS_agn_fsat_*.fits (updates)
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 2  >   EFEDS_AGN_MD10_fsat2_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 4  >   EFEDS_AGN_MD10_fsat4_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 6  >   EFEDS_AGN_MD10_fsat6_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 8  >   EFEDS_AGN_MD10_fsat8_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 12 >  EFEDS_AGN_MD10_fsat12_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 14 >  EFEDS_AGN_MD10_fsat14_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 16 >  EFEDS_AGN_MD10_fsat16_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 20 >  EFEDS_AGN_MD10_fsat20_MODEL.log &
nohup nice -n 19 python 003_0_agn_EFEDS_WTHETA_model.py MD10 10 >  EFEDS_AGN_MD10_fsat10_MODEL.log &

"""
quick turn around verification of the logNlogS
git pull
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_10.fits
python 003_0_agn_EFEDS_WTHETA_logNlogS.py
mv $GIT_AGN_MOCK/figures/MD10/eFEDs/logN_logS_AGN_hard.png $GIT_AGN_MOCK/figures/MD10/eFEDs/Mod1_logN_logS_AGN_hard.png
mv $GIT_AGN_MOCK/figures/MD10/eFEDs/logN_logS_AGN_soft.png $GIT_AGN_MOCK/figures/MD10/eFEDs/Mod1_logN_logS_AGN_soft.png
git add $GIT_AGN_MOCK/figures/MD10/eFEDs/*
git commit -am"new logNlogS"
git push
"""

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

# public link (300G !)
#  cp -r $MD10/cat_AGN_all /data42s/comparat/firefly/mocks/2020-03/MDPL2/

nohup nice -n 19 python run_md_agn.py MD04 all          > run_agn_MD04_all.log & #                      
nohup nice -n 19 python run_md_agn.py MD10 all          > run_agn_MD10_all.log & #                      
nohup nice -n 19 python run_md_agn.py MD40 all          > run_agn_MD40_all.log & #                      
nohup nice -n 19 python run_md_agn.py UNIT_fA1_DIR all  > run_agn_UNIT_fA1_DIR_all.log & #      
nohup nice -n 19 python run_md_agn.py UNIT_fA2_DIR all  > run_agn_UNIT_fA2_DIR_all.log & #      
nohup nice -n 19 python run_md_agn.py UNIT_fA1i_DIR all > run_agn_UNIT_fA1i_DIR_all.log & #  

# create AGN SIMPUT files for eROSITA. Takes 2.5 hours
# Considers a subset of the AGNs in the simulation: FX_soft>2e-17
# script: 003_4_agn_simput.py
# inputs: 768 files in : ${env}/cat_AGN_all/000???.fit 
# outputs: 768 files in : ${env}/cat_AGN_SIMPUT/000???.fit 
# public link (35G)
#  cp -r $MD10/cat_AGN_SIMPUT /data42s/comparat/firefly/mocks/2020-03/MDPL2/
# 
# alternative script, per pixel, more parallel: python 003_4_agn_simput_HPX.py MD10 767
nohup nice -n 19 python 003_4_agn_simput.py MD04          > run_simput_MD04.log &
nohup nice -n 19 python 003_4_agn_simput.py MD10          > run_simput_MD10.log &
nohup nice -n 19 python 003_4_agn_simput.py MD40          > run_simput_MD40.log &
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA1_DIR  > run_simput_UNIT_fA1_DIR.log &
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA2_DIR  > run_simput_UNIT_fA2_DIR.log &
nohup nice -n 19 python 003_4_agn_simput.py UNIT_fA1i_DIR > run_simput_UNIT_fA1i_DIR.log &

# simulate counts with sixte 
# script: $GIT_AGN_MOCK/python/sixte/simulate_agn_only.py
# inputs: 
#    768 catalogues: ${env}/cat_AGN_SIMPUT/000???.fit
#    spectra: ${env}/cat_AGN_SIMPUT/agn_Xspectra/NH??.?_Z?.?_N*.fits
# outputs:
#    /data40s/erosim/eRASS/eRASS8_agn_MD10/???/erass_ccd?_evt.fits
#    /data40s/erosim/eRASS/eRASS8_agn_MD10/???/erass.gti
# concatenated outputs :
#    /data40s/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccd?.fits
#    /data40s/erosim/eRASS/eRASS8_agn_MD10/
# public version
#    cp /data40s/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccd?.fits /data42s/comparat/firefly/mocks/2020-03/MDPL2/AGN_counts/
#    https://firefly.mpe.mpg.de/mocks/2020-03/MDPL2/AGN_counts/

# script for the other simulations are to be created
# execute all commands in this file by copy pasting them into a screen on ds43 of ds54
$GIT_AGN_MOCK/python/sixte/all_sky_simulate_agn_MD10.sh


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
# public link (12G)
#  cp -r $MD10/cat_AGN-MAG_all /data42s/comparat/firefly/mocks/2020-03/MDPL2/

ssh he1s
cd $GIT_AGN_MOCK/python/003_6_log/
pyCONDA
ipython

import os
import glob
import numpy as n
scripts = sorted(n.array(glob.glob("MD04_*.sh")))
for script in scripts:
    os.system('sbatch ' + script)

import os
import glob
import numpy as n
scripts = sorted(n.array(glob.glob("MD10_*.sh")))
for script in scripts:
    os.system('sbatch ' + script)

import os
import glob
import numpy as n
scripts = sorted(n.array(glob.glob("MD40_*.sh")))
for script in scripts:
    os.system('sbatch ' + script)


import os
import glob
import numpy as n
scripts = sorted(n.array(glob.glob("UNIT_fA1i_DIR_*.sh")))
for script in scripts:
    os.system('sbatch ' + script)


import os
import glob
import numpy as n
scripts = sorted(n.array(glob.glob("UNIT_fA2_DIR_*.sh")))
for script in scripts:
    os.system('sbatch ' + script)


import os
import glob
import numpy as n
scripts = sorted(n.array(glob.glob("UNIT_fA1_DIR_*.sh")))
for script in scripts:
    os.system('sbatch ' + script)


# compute SDSS-5 / 4MOST AGN / QSO catalogues for BHM and 4MOST consortium surveys S6 and S8
# script: 003_7_4most_catalog.py
# inputs: 768 files in : ${env}/cat_AGN-MAG_all/000???.fit,
# outputs: 768 files in : ${env}/cat_AGN_4MOST/*_000???.fit 
#   global selection :    |g_lat|>10
#      AGN_IR_000???.fit     ( AGN_type=22 | AGN_type=21 ) < 19
#      AGN_DEEP_000???.fit   g_lon>180 & |ecl_lat|>80 & FX > flux_limit_SNR3 - 0.1
#      AGN_WIDE_000???.fit   g_lon>180 & |ecl_lat|<80 & FX > flux_limit_SNR3 - 0.1
#      LyA_000???.fit        AGN_type==11 & AGN_SDSS_r_magnitude < 23.5 & z > 2.2
#      QSO_000???.fit        AGN_type==11 & AGN_SDSS_r_magnitude < 23.5 & z < 2.2
# 
# public link (12G)
#  cp -r $MD10/cat_AGN-MAG_all /data42s/comparat/firefly/mocks/2020-03/MDPL2/
nohup nice -n 19 python 003_7_4most_catalog.py MD04          > agn_4most_MD04.log &  
nohup nice -n 19 python 003_7_4most_catalog.py MD10          > agn_4most_MD10.log &  
nohup nice -n 19 python 003_7_4most_catalog.py MD40          > agn_4most_MD40.log &  
nohup nice -n 19 python 003_7_4most_catalog.py UNIT_fA1_DIR  > agn_4most_UNIT_fA1_DIR.log &  
nohup nice -n 19 python 003_7_4most_catalog.py UNIT_fA2_DIR  > agn_4most_UNIT_fA2_DIR.log &  
nohup nice -n 19 python 003_7_4most_catalog.py UNIT_fA1i_DIR > agn_4most_UNIT_fA1i_DIR.log &  

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

# MD04  
ls $MD04/cat_AGN-MAG_all/*.fit > md04_4most_agn_all.lis
stilts tcat in=@md04_4most_agn_all.lis ifmt=fits omode=out ofmt=fits out=$MD04/AGN_all.fits
ls $MD04/cat_AGN_4MOST/AGN_DEEP_*.fit > md04_4most_agn_deep.lis
stilts tcat in=@md04_4most_agn_deep.lis ifmt=fits omode=out ofmt=fits out=$MD04/AGN_DEEP_4MOST.fits
ls $MD04/cat_AGN_4MOST/AGN_WIDE_*.fit > md04_4most_agn_wide.lis
stilts tcat in=@md04_4most_agn_wide.lis ifmt=fits omode=out ofmt=fits out=$MD04/AGN_WIDE_4MOST.fits
ls $MD04/cat_AGN_4MOST/AGN_IR_*.fit > md04_4most_agn_ir.lis
stilts tcat in=@md04_4most_agn_ir.lis ifmt=fits omode=out ofmt=fits out=$MD04/AGN_IR_4MOST.fits
ls $MD04/cat_AGN_4MOST/QSO_*.fit > md04_4most_qso.lis
stilts tcat in=@md04_4most_qso.lis ifmt=fits omode=out ofmt=fits out=$MD04/QSO_4MOST.fits
ls $MD04/cat_AGN_4MOST/LyA_*.fit > md04_4most_lya.lis
stilts tcat in=@md04_4most.lis ifmt=fits omode=out ofmt=fits out=$MD04/LyA_4MOST.fits

# MD04  
ls $MD10/cat_AGN-MAG_all/*.fit > 4most_agn_all.lis
stilts tcat in=@4most_agn_all.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_all.fits
ls $MD10/cat_AGN_4MOST/AGN_DEEP_*.fit > 4most_agn_deep.lis
stilts tcat in=@4most_agn_deep.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_DEEP_4MOST.fits
ls $MD10/cat_AGN_4MOST/AGN_WIDE_*.fit > 4most_agn_wide.lis
stilts tcat in=@4most_agn_wide.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_WIDE_4MOST.fits
ls $MD10/cat_AGN_4MOST/AGN_IR_*.fit > 4most_agn_ir.lis
stilts tcat in=@4most_agn_ir.lis ifmt=fits omode=out ofmt=fits out=$MD10/AGN_IR_4MOST.fits
ls $MD10/cat_AGN_4MOST/QSO_*.fit > 4most_agn.lis
stilts tcat in=@4most_agn.lis ifmt=fits omode=out ofmt=fits out=$MD10/QSO_4MOST.fits
ls $MD10/cat_AGN_4MOST/LyA_*.fit > 4most_agn.lis
stilts tcat in=@4most_agn.lis ifmt=fits omode=out ofmt=fits out=$MD10/LyA_4MOST.fits

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
nohup nice -n 19 python 003_7_shared_agn.py MD10 > agn_overlap_matrix_MD10.log & 
nohup nice -n 19 python 003_7_shared_agn.py MD04 > agn_overlap_matrix_MD04.log & 

# Merge the AGN catalogues from MD04 and MD10 into a single 4MOST S6 catalogue
# script: 003_7_merge.py
# inputs: from MD04 and MD10 
#     ${env}/AGN_DEEP_4MOST.fits
#     ${env}/AGN_WIDE_4MOST.fits
#     ${env}/AGN_IR_4MOST.fits
# output: /data42s/comparat/firefly/mocks/2020-03/QMOST/S6_4MOST_R-228-234.fit 
nohup nice -n 19 python 003_7_merge.py > agn_merge_4most_catalogue.log &

# public links 
cp $MD10/AGN_all.fits /data42s/comparat/firefly/mocks/2020-03/MDPL2/MD10_AGN_all.fits
cp $MD04/AGN_all.fits /data42s/comparat/firefly/mocks/2020-03/SMDPL/MD04_AGN_all.fits

# Merge the cosmology QSO catalogue late after adding the footprint bits


