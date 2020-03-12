#####################################################
#
# 4MOST / SDSS-V mock catalogue
#
#####################################################

# creates BCG and CGAL catalogs based on the red galaxy catalogue
# writes catalogues in the 4MOST format
# script: 004_9_4most_catalog.py
# inputs: ${env}/${env}_eRO_CLU_GAL.fit
# outputs: 
#   ${env}/S5_BCG_4MOST.fit
#   ${env}/S5_CGAL_4MOST.fit
nohup nice -n 19 python 004_9_4most_catalog.py MD04 > 004_9_4most_catalog_MD04.log & 
nohup nice -n 19 python 004_9_4most_catalog.py MDa0 > 004_9_4most_catalog_MD10.log & 

# Can only be merges into a full S5 catalogue after the cosmology mocks are made i.e. when the Filament survey mock exists.
# script 004_9_4most_merge.py
# see below in the script
# final file will be /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S5_CLUSTER.fits

######################################
# 4MOST cosmology 
######################################

# preparation scripts 
# tabulate NZ, only run once
# python 005_NZ_tabulateJson.py
# cd $GIT_AGN_MOCK/data/cosmo-4most/
# python 4most_ts_redshifttruth.py --dz=0.1 -n=2 --zmax=0.6


# write QMOST S8 cosmology catalogues and S5 Filament sub survey catalogue
# abundance matching with the NZ tabulated in advance
# script: 005_0_sham_cosmology_catalogs.py
# inputs: 
#    ${env}/cat_GALAXY_all/000???.fit
# outputs: 
#    ${env}/cat_SHAM_COSMO/
#       S5GAL_*
#       BG_*
#       LRG_*
#       ELG_*

# execute the commands in the following script for MD10 and for MD04
sh 005_0_sham_cosmology_catalogs.sh

# Adds K-band magnitude in each file 
# script: 005_1_Kband_magnitudes.py
# inputs: 
#    ${env}/cat_SHAM_COSMO/
#       S5GAL_*
#       BG_*
#       LRG_*
#       ELG_*
# outputs: updates the files
nohup nice -n 19 python 005_1_Kband_magnitudes.py MD10 >  005_1_Kband_mag_MD10.log &
nohup nice -n 19 python 005_1_Kband_magnitudes.py MD04 >  005_1_Kband_mag_MD04.log &

# Adds all magnitudes and the r-band fiber magnitude in each file 
# script: 005_2_all_magnitudes.py
# inputs: 
#    ${env}/cat_SHAM_COSMO/
#       S5GAL_*
#       BG_*
#       LRG_*
#       ELG_*
# outputs: updates the files
nohup nice -n 19 python 005_2_all_magnitudes.py   MD10 >  005_2_all_mag_MD10.log &
nohup nice -n 19 python 005_2_all_magnitudes.py   MD04 >  005_2_all_mag_MD04.log &

# converts to the QMOST format
# script: 005_3_4most_catalog.py
# inputs: 
#    ${env}/cat_SHAM_COSMO/
#       S5GAL_*
#       BG_*
#       LRG_*
#       ELG_*
# outputs: updates the files
nohup nice -n 19 python 005_3_4most_catalog.py    MD10 >  005_3_4most_cat_MD10.log &                                 
nohup nice -n 19 python 005_3_4most_catalog.py    MD04 >  005_3_4most_cat_MD04.log &

###############################
#
# Concatenates ctalogues
#
###############################
ls $MD10/cat_SHAM_COSMO/S5GAL_*.fit > 4most_bgs5.lis
stilts tcat in=@4most_bgs5.lis ifmt=fits omode=out ofmt=fits out=$MD10/FILAMENTS5_4MOST.fits
rm 4most_bgs5.lis
ls $MD10/cat_SHAM_COSMO/BG_*.fit > 4most_bg.lis
stilts tcat in=@4most_bg.lis ifmt=fits omode=out ofmt=fits out=$MD10/BG_4MOST.fits
rm 4most_bg.lis
ls $MD10/cat_SHAM_COSMO/LRG_*.fit > 4most_lrg.lis
stilts tcat in=@4most_lrg.lis ifmt=fits omode=out ofmt=fits out=$MD10/LRG_4MOST.fits
rm 4most_lrg.lis
ls $MD10/cat_SHAM_COSMO/ELG_*.fit > 4most_elg.lis
stilts tcat in=@4most_elg.lis ifmt=fits omode=out ofmt=fits out=$MD10/ELG_4MOST.fits
rm 4most_elg.lis

ls $MD04/cat_SHAM_COSMO/S5GAL_*.fit > 4most_bgs5.lis
stilts tcat in=@4most_bgs5.lis ifmt=fits omode=out ofmt=fits out=$MD04/FILAMENTS5_4MOST.fits
rm 4most_bgs5.lis
ls $MD04/cat_SHAM_COSMO/BG_*.fit > 4most_bg.lis
stilts tcat in=@4most_bg.lis ifmt=fits omode=out ofmt=fits out=$MD04/BG_4MOST.fits
rm 4most_bg.lis
ls $MD04/cat_SHAM_COSMO/LRG_*.fit > 4most_lrg.lis
stilts tcat in=@4most_lrg.lis ifmt=fits omode=out ofmt=fits out=$MD04/LRG_4MOST.fits
rm 4most_lrg.lis
ls $MD04/cat_SHAM_COSMO/ELG_*.fit > 4most_elg.lis
stilts tcat in=@4most_elg.lis ifmt=fits omode=out ofmt=fits out=$MD04/ELG_4MOST.fits
rm 4most_elg.lis

# Add footprint bits
# script: 005_4_add_footprint.py
# inputs: catalogue file (any)    
# outputs: updates the files with the footprint bit

# S5
python 005_4_add_footprint.py MD10 S5_BCG_4MOST.fit
python 005_4_add_footprint.py MD10 S5_CGAL_4MOST.fit
python 005_4_add_footprint.py MD10 FILAMENTS5_4MOST.fits

python 005_4_add_footprint.py MD04 S5_BCG_4MOST.fit
python 005_4_add_footprint.py MD04 S5_CGAL_4MOST.fit
python 005_4_add_footprint.py MD04 FILAMENTS5_4MOST.fits

# S8
python 005_4_add_footprint.py MD10 BG_4MOST.fits
python 005_4_add_footprint.py MD10 LRG_4MOST.fits
python 005_4_add_footprint.py MD10 ELG_4MOST.fits
python 005_4_add_footprint.py MD10 QSO_4MOST.fits
python 005_4_add_footprint.py MD10 LyA_4MOST.fits

python 005_4_add_footprint.py MD04 BG_4MOST.fits
python 005_4_add_footprint.py MD04 LRG_4MOST.fits
python 005_4_add_footprint.py MD04 ELG_4MOST.fits
python 005_4_add_footprint.py MD04 QSO_4MOST.fits
python 005_4_add_footprint.py MD04 LyA_4MOST.fits

# S6
python 005_4_add_footprint.py /data42s/comparat/firefly/mocks/2020-03/QMOST/S6_4MOST_R-228-234.fit

###############################
#
# Merges high-z low-z parts
#
###############################
# S8 cosmology
# public file
# /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_QSO_LyA.fits
# /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_ELG.fits
# /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_BG.fits
# /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_LRG.fits

# S8 Combine MD04 (low Z) and MD10 (high Z)
stilts tpipe in=$MD04/QSO_4MOST.fits ifmt=fits omode=out ofmt=fits out=$MD04/S8_AGN_low_z.fits cmd='select "redshift_estimate<=0.3 "'
stilts tcat ifmt=fits in=$MD10/QSO_4MOST.fits in=$MD10/LyA_4MOST.fits ocmd='select "redshift_estimate>0.3 "' out=$MD10/S8_AGN_hig_z.fit
# summary file
stilts tcat ifmt=fits in=$MD04/S8_AGN_low_z.fits in=$MD10/S8_AGN_hig_z.fit out=/data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_QSO_LyA.fits

stilts tcat ifmt=fits in1=$MD10/BG_4MOST.fits in2=$MD10/BG_4MOST.fits icmd1='select "redshift_estimate>0.3 "' icmd2='select "redshift_estimate<=0.3 "' out=/data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_BG.fits

stilts tcat ifmt=fits in1=$MD10/LRG_4MOST.fits in2=$MD10/LRG_4MOST.fits icmd1='select "redshift_estimate>0.3 "' icmd2='select "redshift_estimate<=0.3 "' out=/data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_LRG.fits

stilts tcat ifmt=fits in1=$MD10/ELG_4MOST.fits in2=$MD10/ELG_4MOST.fits icmd1='select "redshift_estimate>0.3 "' icmd2='select "redshift_estimate<=0.3 "' out=/data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S8_ELG.fits


# Merges BCG and CGAL catalogs into a single 4MOST catalogue
# writes catalogues in the 4MOST format
# script: 004_9_4most_merge.py
# inputs: ${env}/S5_BCG_4MOST.fit, ${env}/S5_CGAL_4MOST.fit, ${env}/FILAMENTS5_4MOST.fit
# outputs: ${env}/S5_4MOST.fit
nohup nice -n 19 python 004_9_4most_merge.py MD04 > 004_9_4most_catalog_MD04.log & 
nohup nice -n 19 python 004_9_4most_merge.py MD10 > 004_9_4most_catalog_MD10.log & 

stilts tpipe in=$MD04/S5_4MOST.fit ifmt=fits omode=out ofmt=fits out=$MD04/S5_4MOST_low_z.fits cmd='select "REDSHIFT_ESTIMATE<0.3"'
stilts tpipe in=$MD10/S5_4MOST.fit ifmt=fits omode=out ofmt=fits out=$MD10/S5_4MOST_hig_z.fit cmd='select "REDSHIFT_ESTIMATE>=0.3"'

stilts tcat ifmt=fits in=$MD04/S5_4MOST_low_z.fits in=$MD10/S5_4MOST_hig_z.fit omode=out out=/data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S5_CLUSTER.fits

# S5 clusters
# public file
# /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S5_CLUSTER.fits
# /data42s/comparat/firefly/mocks/2020-03/QMOST/QMOST_S6_AGN.fits

# compresses the files 
cd /data42s/comparat/firefly/mocks/2020-03/QMOST
gzip -k --rsyncable QMOST_S6_AGN.fits
gzip -k --rsyncable QMOST_S5_CLUSTER.fits
gzip -k --rsyncable QMOST_S8_QSO_LyA.fits
gzip -k --rsyncable QMOST_S8_ELG.fits
gzip -k --rsyncable QMOST_S8_BG.fits
gzip -k --rsyncable QMOST_S8_LRG.fits


### PROVIDE the LIGHT CONE FILES by running

tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.21690.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.21690.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.22170.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.22170.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.22670.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.22670.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.23180.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.23180.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.23690.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.23690.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.24230.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.24230.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.24770.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.24770.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.25320.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.25320.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.25890.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.25890.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.26470.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.26470.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.27060.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.27060.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.27670.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.27670.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.28290.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.28290.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.28920.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.28920.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.29570.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.29570.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.30230.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.30230.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.30910.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.30910.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.31600.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.31600.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.32310.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.32310.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.33030.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.33030.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.33770.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.33770.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.34530.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.34530.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.35300.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.35300.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.36090.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.36090.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.36900.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.36900.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.37730.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.37730.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.38570.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.38570.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.39440.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.39440.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.40320.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.40320.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.41230.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.41230.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.42150.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.42150.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.43090.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.43090.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.44060.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.44060.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.45050.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.45050.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.46050.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.46050.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.47090.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.47090.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.48140.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.48140.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.49220.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.49220.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.50320.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.50320.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.51450.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.51450.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.52600.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.52600.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.53780.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.53780.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.54980.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.54980.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.56220.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.56220.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.57470.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.57470.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.58760.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.58760.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.60080.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.60080.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.61420.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.61420.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.62800.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.62800.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.64210.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.64210.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.65650.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.65650.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.67120.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.67120.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.68620.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.68620.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.70160.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.70160.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.71730.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.71730.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.73330.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.73330.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.74980.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.74980.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.76660.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.76660.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.78370.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.78370.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.80130.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.80130.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.81920.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.81920.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.83760.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.83760.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.85640.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.85640.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.87550.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.87550.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.89510.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.89510.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.91520.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.91520.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.93570.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.93570.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.95670.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.95670.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_0.97810.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_0.97810.fits
tar -czf /home/comparat/data2/firefly/mocks/2020-03/MDPL2/LC/all_1.00000.fits.tar.gz /data37s/simulation_1/MD/MD_1.0Gpc/fits/all_1.00000.fits


