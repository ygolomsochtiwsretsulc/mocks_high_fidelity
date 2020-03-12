#############################################
#############################################
#############################################
#
# CLUSTER model
#
#############################################
#############################################
#############################################

# computes cluster properties up to z=0.4 for MD04 and 1.2 (aexp=0.45) for the other simulations
# writes SIMPUT images (option yes, no)
# script: 004_0_cluster.py
# inputs: ${env}/fits/all_?.?????.fits, ${env}/fits/all_?.?????_coordinates.fits, ${env}/fits/all_?.?????_galaxy.fits
# outputs: ${env}/fits/all_?.?????_CLU.fits

# Serial way to run it :
nohup nice -n 19  python run_md_cluster.py MD04 no >  cluster_run_MD04_no.log &                    
nohup nice -n 19  python run_md_cluster.py MD10 no >  cluster_run_MD10_no.log &                    
nohup nice -n 19  python run_md_cluster.py MD40 no >  cluster_run_MD40_no.log &                    
nohup nice -n 19  python run_md_cluster.py UNIT_fA1_DIR  no >  cluster_run_UNIT_fA1_DIR_no.log &  
nohup nice -n 19  python run_md_cluster.py UNIT_fA2_DIR  no >  cluster_run_UNIT_fA2_DIR_no.log &    
nohup nice -n 19  python run_md_cluster.py UNIT_fA1i_DIR no >  cluster_run_UNIT_fA1i_DIR_no.log &  

nohup nice -n 19  python run_md_cluster.py MD04 yes >  cluster_run_MD04_yes.log &                   
nohup nice -n 19  python run_md_cluster.py MD10 yes >  cluster_run_MD10_yes.log &                   
nohup nice -n 19  python run_md_cluster.py MD40 yes >  cluster_run_MD40.log &                    
nohup nice -n 19  python run_md_cluster.py UNIT_fA1_DIR  yes >  cluster_run_UNIT_fA1_DIR.log &  
nohup nice -n 19  python run_md_cluster.py UNIT_fA2_DIR  yes >  cluster_run_UNIT_fA2_DIR.log &    
nohup nice -n 19  python run_md_cluster.py UNIT_fA1i_DIR yes >  cluster_run_UNIT_fA1i_DIR.log &  

#####################################################
#
# File by file (parallel way to run it)
#
#####################################################

############## MD10 ################
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.45050 yes >  cluster_run_MD10_0.45050_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.46050 yes >  cluster_run_MD10_0.46050_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.47090 yes >  cluster_run_MD10_0.47090_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.48140 yes >  cluster_run_MD10_0.48140_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.49220 yes >  cluster_run_MD10_0.49220_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.50320 yes >  cluster_run_MD10_0.50320_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.51450 yes >  cluster_run_MD10_0.51450_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.52600 yes >  cluster_run_MD10_0.52600_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.53780 yes >  cluster_run_MD10_0.53780_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.54980 yes >  cluster_run_MD10_0.54980_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.56220 yes >  cluster_run_MD10_0.56220_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.57470 yes >  cluster_run_MD10_0.57470_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.58760 yes >  cluster_run_MD10_0.58760_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.60080 yes >  cluster_run_MD10_0.60080_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.61420 yes >  cluster_run_MD10_0.61420_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.62800 yes >  cluster_run_MD10_0.62800_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.64210 yes >  cluster_run_MD10_0.64210_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.65650 yes >  cluster_run_MD10_0.65650_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.67120 yes >  cluster_run_MD10_0.67120_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.68620 yes >  cluster_run_MD10_0.68620_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.70160 yes >  cluster_run_MD10_0.70160_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.71730 yes >  cluster_run_MD10_0.71730_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.73330 yes >  cluster_run_MD10_0.73330_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.74980 yes >  cluster_run_MD10_0.74980_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.76660 yes >  cluster_run_MD10_0.76660_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.78370 yes >  cluster_run_MD10_0.78370_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.80130 yes >  cluster_run_MD10_0.80130_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.81920 yes >  cluster_run_MD10_0.81920_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.83760 yes >  cluster_run_MD10_0.83760_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.85640 yes >  cluster_run_MD10_0.85640_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.87550 yes >  cluster_run_MD10_0.87550_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.89510 yes >  cluster_run_MD10_0.89510_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.91520 yes >  cluster_run_MD10_0.91520_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.93570 yes >  cluster_run_MD10_0.93570_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.95670 yes >  cluster_run_MD10_0.95670_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_0.97810 yes >  cluster_run_MD10_0.97810_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD10 all_1.00000 yes >  cluster_run_MD10_1.00000_yes.log &

############## MD04 ################
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.45050 yes >  cluster_run_MD40_0.45050_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.46050 yes >  cluster_run_MD40_0.46050_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.47090 yes >  cluster_run_MD40_0.47090_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.48140 yes >  cluster_run_MD40_0.48140_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.49220 yes >  cluster_run_MD40_0.49220_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.50320 yes >  cluster_run_MD40_0.50320_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.51450 yes >  cluster_run_MD40_0.51450_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.52600 yes >  cluster_run_MD40_0.52600_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.53780 yes >  cluster_run_MD40_0.53780_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.54980 yes >  cluster_run_MD40_0.54980_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.56220 yes >  cluster_run_MD40_0.56220_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.57470 yes >  cluster_run_MD40_0.57470_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.58760 yes >  cluster_run_MD40_0.58760_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.60080 yes >  cluster_run_MD40_0.60080_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.61420 yes >  cluster_run_MD40_0.61420_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.62800 yes >  cluster_run_MD40_0.62800_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.64210 yes >  cluster_run_MD40_0.64210_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.65650 yes >  cluster_run_MD40_0.65650_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.67120 yes >  cluster_run_MD40_0.67120_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.68620 yes >  cluster_run_MD40_0.68620_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.70160 yes >  cluster_run_MD40_0.70160_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.71730 yes >  cluster_run_MD40_0.71730_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.73330 yes >  cluster_run_MD40_0.73330_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.74980 yes >  cluster_run_MD40_0.74980_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.76660 yes >  cluster_run_MD40_0.76660_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.78370 yes >  cluster_run_MD40_0.78370_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.80130 yes >  cluster_run_MD40_0.80130_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.81920 yes >  cluster_run_MD40_0.81920_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.83760 yes >  cluster_run_MD40_0.83760_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.85640 yes >  cluster_run_MD40_0.85640_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.87550 yes >  cluster_run_MD40_0.87550_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.89510 yes >  cluster_run_MD40_0.89510_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.91520 yes >  cluster_run_MD40_0.91520_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.93570 yes >  cluster_run_MD40_0.93570_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.95670 yes >  cluster_run_MD40_0.95670_yes.log &      
nohup nice -n 19 python3 004_0_cluster.py MD40 all_0.97810 yes >  cluster_run_MD40_0.97810_yes.log &      

############## MD40 ################
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.70030 yes >  cluster_run_MD04_0.70030_yes.log &    
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.70640 yes >  cluster_run_MD04_0.70640_yes.log &     
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.71240 yes >  cluster_run_MD04_0.71240_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.71830 yes >  cluster_run_MD04_0.71830_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.72430 yes >  cluster_run_MD04_0.72430_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.73630 yes >  cluster_run_MD04_0.73630_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.74230 yes >  cluster_run_MD04_0.74230_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.75440 yes >  cluster_run_MD04_0.75440_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.76030 yes >  cluster_run_MD04_0.76030_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.77240 yes >  cluster_run_MD04_0.77240_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.77830 yes >  cluster_run_MD04_0.77830_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.78730 yes >  cluster_run_MD04_0.78730_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.79040 yes >  cluster_run_MD04_0.79040_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.80230 yes >  cluster_run_MD04_0.80230_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.80840 yes >  cluster_run_MD04_0.80840_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.81440 yes >  cluster_run_MD04_0.81440_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.81730 yes >  cluster_run_MD04_0.81730_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.82340 yes >  cluster_run_MD04_0.82340_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.82630 yes >  cluster_run_MD04_0.82630_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.83240 yes >  cluster_run_MD04_0.83240_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.83530 yes >  cluster_run_MD04_0.83530_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.84140 yes >  cluster_run_MD04_0.84140_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.84430 yes >  cluster_run_MD04_0.84430_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.85040 yes >  cluster_run_MD04_0.85040_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.85330 yes >  cluster_run_MD04_0.85330_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.85940 yes >  cluster_run_MD04_0.85940_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.86230 yes >  cluster_run_MD04_0.86230_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.86840 yes >  cluster_run_MD04_0.86840_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.87130 yes >  cluster_run_MD04_0.87130_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.87730 yes >  cluster_run_MD04_0.87730_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.88030 yes >  cluster_run_MD04_0.88030_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.88630 yes >  cluster_run_MD04_0.88630_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.88930 yes >  cluster_run_MD04_0.88930_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.89530 yes >  cluster_run_MD04_0.89530_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.89840 yes >  cluster_run_MD04_0.89840_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.90430 yes >  cluster_run_MD04_0.90430_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.90740 yes >  cluster_run_MD04_0.90740_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.92515 yes >  cluster_run_MD04_0.92515_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.95600 yes >  cluster_run_MD04_0.95600_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_0.97071 yes >  cluster_run_MD04_0.97071_yes.log &
nohup nice -n 19 python3 004_0_cluster.py MD04 all_1.00000 yes >  cluster_run_MD04_1.00000_yes.log &

#####################################################
#
# Merge the cluster catalogue
#
#####################################################

ls $MD04/fits/*_CLU.fits > fit_list_eRo_CLUMD04.list
stilts tcat in=@fit_list_eRo_CLUMD04.list ifmt=fits omode=out ofmt=fits out=$MD04/MD04_eRO_CLU.fit

ls $MD40/fits/*_CLU.fits > fit_list_eRo_CLUMD40.list
stilts tcat in=@fit_list_eRo_CLUMD40.list ifmt=fits omode=out ofmt=fits out=$MD40/MD40_eRO_CLU.fit

ls $MD10/fits/*_CLU.fits > fit_list_eRo_CLUMD10.list
stilts tcat in=@fit_list_eRo_CLUMD10.list ifmt=fits omode=out ofmt=fits out=$MD10/MD10_eRO_CLU.fit

###########################################
#
# FIGURES of clusters
#
###########################################

# plots the mass (M500c) vs redshift for the simulation and the data used by the cluster profile generator
# script: 004_5_plot_clusters_M500-z.py
# inputs: ${env}/${env}_eRO_CLU.fit
# outputs: 
#   $GIT_AGN_MOCK/figures/${env}/clustersM500c-z-1em14.png (for a given flux limit)
#   $GIT_AGN_MOCK/figures/${env}/clustersM500c-zpng (full simulation)
cd $GIT_AGN_MOCK/python/
nohup nice -n 19  python 004_5_plot_clusters_M500-z.py MD10 200c  >  plot_004_5_M500z_MD04.log &
nohup nice -n 19  python 004_5_plot_clusters_M500-z.py MD04 200c  >  plot_004_5_M500z_MD10.log &

# plots the logNlogS
# script: plot_clusters_logNlogS.py
# inputs: ${env}/${env}_eRO_CLU.fit
# outputs: 
# figures per simulation in $GIT_AGN_MOCK/figures/${env}/clusters/
#   logNlogS_2e13.png logNlogS_2e13.txt
#   logNlogS_5e13.png logNlogS_5e13.txt
#   logNlogS_1e14.png logNlogS_1e14.txt
#  summary figures (can accomodate all simulations, but currently only MD10)
#   $GIT_AGN_MOCK/figures/logNlogS_cluster_2e13.png
#   $GIT_AGN_MOCK/figures/logNlogS_cluster_5e13.png
#   $GIT_AGN_MOCK/figures/logNlogS_cluster_1e14.png
cd $GIT_AGN_MOCK/python/figures/
nohup nice -n 19 python plot_clusters_logNlogS.py MD10                        > lognlogs_MD10.log           &
nohup nice -n 19 python plot_clusters_logNlogS.py MD04                        > lognlogs_MD04.log           &
nohup nice -n 19 python plot_clusters_logNlogS.py MD40                        > lognlogs_MD40.log           &
nohup nice -n 19 python plot_clusters_logNlogS.py UNIT_fA1i_DIR               > lognlogs_UNIT_fA1i_DIR.log  &
nohup nice -n 19 python plot_clusters_logNlogS.py UNIT_fA1_DIR                > lognlogs_UNIT_fA1_DIR.log   &
nohup nice -n 19 python plot_clusters_logNlogS.py UNIT_fA2_DIR                > lognlogs_UNIT_fA2_DIR.log   &

# plots the X-ray luminosity function
# script: plot_clusters_XLF.py
# inputs: ${env}/${env}_eRO_CLU.fit
# outputs: 
# figure in $GIT_AGN_MOCK/figures/${env}/clusters/XLF.png
cd $GIT_AGN_MOCK/python/figures/
nohup nice -n 19 python plot_clusters_XLF.py MD10                             > xlf_MD10.log           &
nohup nice -n 19 python plot_clusters_XLF.py MD04                             > xlf_MD04.log           &
nohup nice -n 19 python plot_clusters_XLF.py MD40                             > xlf_MD40.log           &
nohup nice -n 19 python plot_clusters_XLF.py UNIT_fA1i_DIR                    > xlf_UNIT_fA1i_DIR.log  &
nohup nice -n 19 python plot_clusters_XLF.py UNIT_fA1_DIR                     > xlf_UNIT_fA1_DIR.log   &
nohup nice -n 19 python plot_clusters_XLF.py UNIT_fA2_DIR                     > xlf_UNIT_fA2_DIR.log   &

# plots the scaling relations LX-kT and M-kT
# script: plot_clusters_scaling_relations.py
# inputs: ${env}/${env}_eRO_CLU.fit
# outputs: 
# figure in $GIT_AGN_MOCK/figures/${env}/clusters/
#   LX-KT-Z.png
#   M500-KT-Z.png
cd $GIT_AGN_MOCK/python/figures/
nohup nice -n 19 python plot_clusters_scaling_relations.py MD04               > scaling_relations_MD04.log           &
nohup nice -n 19 python plot_clusters_scaling_relations.py MD10               > scaling_relations_MD10.log           &
nohup nice -n 19 python plot_clusters_scaling_relations.py MD40               > scaling_relations_MD40.log           &
nohup nice -n 19 python plot_clusters_scaling_relations.py UNIT_fA1i_DIR      > scaling_relations_UNIT_fA1i_DIR.log  &
nohup nice -n 19 python plot_clusters_scaling_relations.py UNIT_fA1_DIR       > scaling_relations_UNIT_fA1_DIR.log   &
nohup nice -n 19 python plot_clusters_scaling_relations.py UNIT_fA2_DIR       > scaling_relations_UNIT_fA2_DIR.log   &

# creates simput catalogues and cluster catalogues with a flux limit
# script: 004_6_clusters_simput.py
# inputs: ${env}/fits/all_?.?????_CLU.fits
# outputs: ${env}/cat_CLU_SIMPUT/c_000???_N_?.fit, ${env}/cat_eRO_CLU/000???.fit
nohup python 004_6_clusters_simput.py MD04          > run_CLU_simput_MD04.log & #                       
nohup python 004_6_clusters_simput.py MD10          > run_CLU_simput_MD10.log & #                       
nohup python 004_6_clusters_simput.py MD40          > run_CLU_simput_MD40.log & #                       

# ONLY RUN ONCE
# computes cluster X-ray spectra
# ds54 or ds43 to have pyxspec
# script: 004_7_clusters_create_xpectra.py
# inputs: env, pyXspec
# outputs: ${env}/cat_CLU_SIMPUT/cluster_Xspectra/cluster_spectrum_10kT_????_100z_????.fits
# nohup nice -n 19 python 004_7_clusters_create_xpectra.py MD10 > run_CLU_7_MD10.log &
# nohup nice -n 19 python 004_7_clusters_create_xpectra.py MD04 > run_CLU_7_MD04.log &

cd $GIT_AGN_MOCK/python/sixte/

# simulate counts with sixte 
# script: $GIT_AGN_MOCK/python/sixte/simulate_cluster_only.py
# inputs: 
#    768 catalogues: ${env}/cat_CLU_SIMPUT/000???.fit
#    spectra: ${env}/cat_CLU_SIMPUT/cluster_Xspectra/NH??.?_Z?.?_N*.fits
#    images: ${env}/cat_CLU_SIMPUT/cluster_images/*
# outputs:
#    /data40s/erosim/eRASS/eRASS8_cluster_${env}/???/erass_ccd?_evt.fits
#    /data40s/erosim/eRASS/eRASS8_cluster_${env}/???/erass.gti
# concatenated outputs :
#    /data40s/erosim/eRASS/eRASS8_cluster_${env}/simulated_photons_ccd?.fits
#    /data40s/erosim/eRASS/eRASS8_cluster_${env}/
# public version
#    cp /data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd?.fits /data42s/comparat/firefly/mocks/2020-03/MDPL2/CLUSTER_counts/
#    cp /data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd?.fits /data42s/comparat/firefly/mocks/2020-03/HMD/CLUSTER_counts/
#    cp /data40s/erosim/eRASS/eRASS8_cluster_MD04/simulated_photons_ccd?.fits /data42s/comparat/firefly/mocks/2020-03/SMDPL/CLUSTER_counts/
#    https://firefly.mpe.mpg.de/mocks/2020-03/MDPL2/CLUSTER_counts/

# execute all commands in 
# execute all commands in this file by copy pasting them into a screen on ds43 of ds54
$GIT_AGN_MOCK/python/sixte/all_sky_simulate_cluster_MD10.sh
$GIT_AGN_MOCK/python/sixte/all_sky_simulate_cluster_MD40.sh
$GIT_AGN_MOCK/python/sixte/all_sky_simulate_cluster_MD04.sh



#####################################################
#
# Galaxies in clusters
#
#####################################################

# ONLY ON DS43 OR 54
# due to healpy dependency

# STEP 1 
# runs all the steps below in a series for each shell
# 
# script: 004_1_galaxies_around_clusters.py
# inputs: env, basename, delta_crit
#    e.g. 'MD04' "all_0.62840" '200c'
# outputs: 
#       runs the scripts below in order

# STEP 2
# identifies galaxies within a certain distance from a halo in the cluster file
# script: 004_2_cluster_galaxies.py
# inputs: ${env}/fits/all_?.?????_CLU.fits, ${env}/fits/all_?.?????.fits, ${env}/fits/all_?.?????_coordinates.fits, ${env}/fits/all_?.?????_galaxy.fits
# outputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits

# STEP 3
# identifies the quiescent galaxies around the clusters
# corrects the star formation rate values
# script: 004_3_cluster_red_galaxies.py
# inputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits
# outputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits (add more columns)

# STEP 4
# Abundance matching to the Ricci et al. 2018 luminosity function fo obtain rest frame r-band magnitudes
# redmapper red-sequence model to give colors to quiescent galaxies
# estimates richness (still in dev.)
# script: 004_4_red_sequence.py
# inputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits
# outputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits (add more columns)

#####################################################
#
# RUN file by file (parallel way to run it)
#
#####################################################

########### MD04 ###########
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.70030 200c > cluster_004_galaxies_MD04_0.70030.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.70640 200c > cluster_004_galaxies_MD04_0.70640.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.71240 200c > cluster_004_galaxies_MD04_0.71240.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.71830 200c > cluster_004_galaxies_MD04_0.71830.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.72430 200c > cluster_004_galaxies_MD04_0.72430.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.73630 200c > cluster_004_galaxies_MD04_0.73630.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.74230 200c > cluster_004_galaxies_MD04_0.74230.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.75440 200c > cluster_004_galaxies_MD04_0.75440.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.76030 200c > cluster_004_galaxies_MD04_0.76030.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.77240 200c > cluster_004_galaxies_MD04_0.77240.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.77830 200c > cluster_004_galaxies_MD04_0.77830.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.78730 200c > cluster_004_galaxies_MD04_0.78730.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.79040 200c > cluster_004_galaxies_MD04_0.79040.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.80230 200c > cluster_004_galaxies_MD04_0.80230.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.80840 200c > cluster_004_galaxies_MD04_0.80840.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.81440 200c > cluster_004_galaxies_MD04_0.81440.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.81730 200c > cluster_004_galaxies_MD04_0.81730.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.82340 200c > cluster_004_galaxies_MD04_0.82340.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.82630 200c > cluster_004_galaxies_MD04_0.82630.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.83240 200c > cluster_004_galaxies_MD04_0.83240.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.83530 200c > cluster_004_galaxies_MD04_0.83530.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.84140 200c > cluster_004_galaxies_MD04_0.84140.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.84430 200c > cluster_004_galaxies_MD04_0.84430.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.85040 200c > cluster_004_galaxies_MD04_0.85040.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.85330 200c > cluster_004_galaxies_MD04_0.85330.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.85940 200c > cluster_004_galaxies_MD04_0.85940.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.86230 200c > cluster_004_galaxies_MD04_0.86230.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.86840 200c > cluster_004_galaxies_MD04_0.86840.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.87130 200c > cluster_004_galaxies_MD04_0.87130.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.87730 200c > cluster_004_galaxies_MD04_0.87730.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.88030 200c > cluster_004_galaxies_MD04_0.88030.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.88630 200c > cluster_004_galaxies_MD04_0.88630.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.88930 200c > cluster_004_galaxies_MD04_0.88930.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.89530 200c > cluster_004_galaxies_MD04_0.89530.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.89840 200c > cluster_004_galaxies_MD04_0.89840.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.90430 200c > cluster_004_galaxies_MD04_0.90430.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.90740 200c > cluster_004_galaxies_MD04_0.90740.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.92515 200c > cluster_004_galaxies_MD04_0.92515.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.95600 200c > cluster_004_galaxies_MD04_0.95600.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_0.97071 200c > cluster_004_galaxies_MD04_0.97071.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD04 all_1.00000 200c > cluster_004_galaxies_MD04_1.00000.log &

########### MD10 ###########
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.45050 200c > cluster_004_gal_MD10_0.45050.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.46050 200c > cluster_004_gal_MD10_0.46050.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.47090 200c > cluster_004_gal_MD10_0.47090.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.48140 200c > cluster_004_gal_MD10_0.48140.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.49220 200c > cluster_004_gal_MD10_0.49220.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.50320 200c > cluster_004_gal_MD10_0.50320.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.51450 200c > cluster_004_gal_MD10_0.51450.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.52600 200c > cluster_004_gal_MD10_0.52600.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.53780 200c > cluster_004_gal_MD10_0.53780.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.54980 200c > cluster_004_gal_MD10_0.54980.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.56220 200c > cluster_004_gal_MD10_0.56220.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.57470 200c > cluster_004_gal_MD10_0.57470.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.58760 200c > cluster_004_gal_MD10_0.58760.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.60080 200c > cluster_004_gal_MD10_0.60080.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.61420 200c > cluster_004_gal_MD10_0.61420.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.62800 200c > cluster_004_gal_MD10_0.62800.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.64210 200c > cluster_004_gal_MD10_0.64210.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.65650 200c > cluster_004_gal_MD10_0.65650.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.67120 200c > cluster_004_gal_MD10_0.67120.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.68620 200c > cluster_004_gal_MD10_0.68620.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.70160 200c > cluster_004_gal_MD10_0.70160.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.71730 200c > cluster_004_gal_MD10_0.71730.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.73330 200c > cluster_004_gal_MD10_0.73330.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.74980 200c > cluster_004_gal_MD10_0.74980.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.76660 200c > cluster_004_gal_MD10_0.76660.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.78370 200c > cluster_004_gal_MD10_0.78370.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.80130 200c > cluster_004_gal_MD10_0.80130.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.81920 200c > cluster_004_gal_MD10_0.81920.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.83760 200c > cluster_004_gal_MD10_0.83760.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.85640 200c > cluster_004_gal_MD10_0.85640.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.87550 200c > cluster_004_gal_MD10_0.87550.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.89510 200c > cluster_004_gal_MD10_0.89510.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.91520 200c > cluster_004_gal_MD10_0.91520.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.93570 200c > cluster_004_gal_MD10_0.93570.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.95670 200c > cluster_004_gal_MD10_0.95670.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_0.97810 200c > cluster_004_gal_MD10_0.97810.log &
nohup nice -n 19 python 004_1_galaxies_around_clusters.py MD10 all_1.00000 200c > cluster_004_gal_MD10_1.00000.log &

###########################################
#
# FIGURES of galaxies in clusters
#
###########################################

# plots the luminosity function for each redshift shell
# script: 004_5_plot_clusters.py
# inputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits
# outputs: $GIT_AGN_MOCK/figures/${env}/clusters/galaxies/cluster_galaxy_LF_z_*.png
cd $GIT_AGN_MOCK/python/
nohup nice -n 19  python 004_5_plot_clusters.py MD04 200c  >  plot_004_5_MD04.log &
nohup nice -n 19  python 004_5_plot_clusters.py MD10 200c  >  plot_004_5_MD10.log &

# plots the distribution of colors: red sequence model and its realization in each redshift shell 
# script: 004_5_plot_clusters_colors.py
# inputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits
# outputs: $GIT_AGN_MOCK/figures/${env}/clusters/galaxies/red_sequence_*.png
cd $GIT_AGN_MOCK/python/
nohup nice -n 19  python 004_5_plot_clusters_colors.py MD04 >  plot_004_5_color_MD04.log &
nohup nice -n 19  python 004_5_plot_clusters_colors.py MD10 >  plot_004_5_color_MD10.log &

# plots the richness vs. mass (M200c) in each redshift shell and writes the data points for fitting
# script: 004_5_plot_clusters_colors.py
# inputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits
# outputs: 
#   $GIT_AGN_MOCK/figures/${env}/clusters/galaxies/richness_mass_*.png
#   $GIT_AGN_MOCK/figures/${env}/clusters/galaxies/richness_mass_*.txt 
cd $GIT_AGN_MOCK/python/
nohup nice -n 19  python 004_5_plot_clusters_richness.py MD10 200c  >  plot_004_5_richness_MD04.log &
nohup nice -n 19  python 004_5_plot_clusters_richness.py MD04 200c  >  plot_004_5_richness_MD10.log &

# tries to fit of the richness scaling relation :
# plots the richness vs. mass (M200c) in each redshift shell and writes the data points for fitting
# script: 004_5_plot_clusters_colors.py
# inputs:  $GIT_AGN_MOCK/figures/${env}/clusters/galaxies/richness_mass_*.txt
# outputs: 
#   set of parameters
#   $GIT_AGN_MOCK/figures/${env}/clusters/galaxies/richness_mass_evolution_1e14.txt and richness_mass_evolution_1e14.png
cd $GIT_AGN_MOCK/python/
nohup nice -n 19  python 004_5_plot_clusters_richness_fits.py MD10 200c  >  plot_004_5_richness_MD04_fits.log &
nohup nice -n 19  python 004_5_plot_clusters_richness_fits.py MD04 200c  >  plot_004_5_richness_MD10_fits.log &

# In progress ...
# Tries to compares galaxy LF with SPIDERS
# script: 004_5_plot_clusters_w_SPIDERS.py
# inputs: ${env}/fits/all_?.?????_galaxiesAroundClusters.fits
# outputs: 
#   $GIT_AGN_MOCK/figures/${env}/galaxies/spiders_cluster_galaxy_LF_z_*.png 
nohup nice -n 19  python 004_5_plot_clusters_w_SPIDERS.py MD04 200c  >  plot_004_5_comparison_spiders_MD04.log &
nohup nice -n 19  python 004_5_plot_clusters_w_SPIDERS.py MD10 200c  >  plot_004_5_comparison_spiders_MD10.log &


#####################################################
#
# Merge the catalogue of galaxies in clusters
# r-band magnitude observed brigther than 24.5
#
#####################################################
 
ls $MD10/fits/*_galaxiesAroundClusters.fit > fit_list_galaxiesAroundClusters_md10.list
stilts tcat in=@fit_list_galaxiesAroundClusters_md10.list ifmt=fits icmd='select "mag_r<24.5"' omode=out ofmt=fits out=$MD10/MD10_eRO_CLU_GAL.fit

ls $MD04/fits/*_galaxiesAroundClusters.fit > fit_list_galaxiesAroundClusters.list
stilts tcat in=@fit_list_galaxiesAroundClusters.list ifmt=fits icmd='select "mag_r<24.5"' omode=out ofmt=fits out=$MD04/MD04_eRO_CLU_GAL.fit

# plots the fraction of red galaxies vs radius
# script: plot_red_fraction_vs_cluster_radius.py
# inputs: ${env}/${env}_eRO_CLU.fit,  ${env}/${env}_eRO_CLU_GAL.fit
# outputs: 
#   $GIT_AGN_MOCK/figures/${env}/frac_QU_radius.png 
cd $GIT_AGN_MOCK/python/figures/
nohup nice -n 19  python plot_red_fraction_vs_cluster_radius.py MD04 >  plot_clusters_MD04.log &                    
nohup nice -n 19  python plot_red_fraction_vs_cluster_radius.py MD10 >  plot_clusters_MD10.log &                    

