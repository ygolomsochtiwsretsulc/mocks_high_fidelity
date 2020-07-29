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
cd $GIT_AGN_MOCK/python/

nohup python run_md_cluster.py MD10 yes 20. 1.0 0 2.0 >    logs/cluster_run_MD10_yes_20_10_0.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.9 0 2.0 >    logs/cluster_run_MD10_yes_20_09_0.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.8 0 2.0 >    logs/cluster_run_MD10_yes_20_08_0.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.7 0 2.0 >    logs/cluster_run_MD10_yes_20_07_0.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.6 0 2.0 >    logs/cluster_run_MD10_yes_20_06_0.log &                    

nohup python run_md_cluster.py MD10 yes 20. 0.8 m20 2.0 >  logs/cluster_run_MD10_yes_20_08_m20.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.8 p20 2.0 >  logs/cluster_run_MD10_yes_20_08_p20.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.8 m40 2.0 >  logs/cluster_run_MD10_yes_20_08_m40.log &                    
nohup python run_md_cluster.py MD10 yes 20. 0.8 p40 2.0 >  logs/cluster_run_MD10_yes_20_08_p40.log &                    

nohup python run_md_cluster.py MD10 yes 2.  0.9 0 4.0 >    logs/cluster_run_MD10_yes_2_09_0.log &                    
nohup python run_md_cluster.py MD10 yes 2.  0.8 0 4.0 >    logs/cluster_run_MD10_yes_2_08_0.log &                    
nohup python run_md_cluster.py MD10 yes 2.  0.7 0 4.0 >    logs/cluster_run_MD10_yes_2_07_0.log &                    

nohup python run_md_cluster.py MD40 yes 2.  0.9 0 4.0 >    logs/cluster_run_MD40_yes_2_09_0.log &                    
nohup python run_md_cluster.py MD40 yes 2.  0.8 0 4.0 >    logs/cluster_run_MD40_yes_2_08_0.log &                    
nohup python run_md_cluster.py MD40 yes 2.  0.7 0 4.0 >    logs/cluster_run_MD40_yes_2_07_0.log &                    

nohup python run_md_cluster.py MD40 yes 20.  0.9 0 4.0 >   logs/cluster_run_MD40_yes_20_09_0.log &                    
nohup python run_md_cluster.py MD40 yes 20.  0.8 0 4.0 >   logs/cluster_run_MD40_yes_20_08_0.log &                    
nohup python run_md_cluster.py MD40 yes 20.  0.7 0 4.0 >   logs/cluster_run_MD40_yes_20_07_0.log &                    

nohup python run_md_cluster.py MD04 yes 20.  0.9 0 0.4 >   logs/cluster_run_MD04_yes_2_09_0.log &                  
nohup python run_md_cluster.py MD04 yes 20.  0.8 0 0.4 >   logs/cluster_run_MD04_yes_2_08_0.log &                  
nohup python run_md_cluster.py MD04 yes 20.  0.7 0 0.4 >   logs/cluster_run_MD04_yes_2_07_0.log &                  

#####################################################
#
# Merge the cluster catalogue
#
#####################################################

cd $GIT_AGN_MOCK/python/

ls $MD04/fits/all_?.?????_CLU_b7_CM_0_pixS_20.0.fits > lists/fit_list_eRO_MD04_b7_CM_0_pixS_20.0.fits
ls $MD04/fits/all_?.?????_CLU_b8_CM_0_pixS_20.0.fits > lists/fit_list_eRO_MD04_b8_CM_0_pixS_20.0.fits
ls $MD04/fits/all_?.?????_CLU_b9_CM_0_pixS_20.0.fits > lists/fit_list_eRO_MD04_b9_CM_0_pixS_20.0.fits

stilts tcat in=@lists/fit_list_eRO_MD04_b7_CM_0_pixS_20.0.fits ifmt=fits omode=out ofmt=fits out=$MD04/MD04_eRO_CLU_b7_CM_0_pixS_20.0.fits
stilts tcat in=@lists/fit_list_eRO_MD04_b8_CM_0_pixS_20.0.fits ifmt=fits omode=out ofmt=fits out=$MD04/MD04_eRO_CLU_b8_CM_0_pixS_20.0.fits
stilts tcat in=@lists/fit_list_eRO_MD04_b9_CM_0_pixS_20.0.fits ifmt=fits omode=out ofmt=fits out=$MD04/MD04_eRO_CLU_b9_CM_0_pixS_20.0.fits


#####################################################
#
# SIMPUT files and sixte run
#
#####################################################

cd $GIT_AGN_MOCK/python/
# creates simput catalogues and cluster catalogues with a flux limit
# script: 004_6_clusters_simput.py
# inputs: ${env}/fits/all_?.?????_CLU.fits
# outputs: ${env}/cat_CLU_SIMPUT/c_000???_N_?.fit, ${env}/cat_eRO_CLU/000???.fit

nohup python 004_6_clusters_simput.py MD10 20. 1.0 0 >    logs/cluster_run_MD10_20_10_0.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.9 0 >    logs/cluster_run_MD10_20_09_0.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.8 0 >    logs/cluster_run_MD10_20_08_0.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.7 0 >    logs/cluster_run_MD10_20_07_0.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.6 0 >    logs/cluster_run_MD10_20_06_0.log &                    

nohup python 004_6_clusters_simput.py MD10 20. 0.8 m20 >  logs/cluster_run_MD10_20_08_m20.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.8 p20 >  logs/cluster_run_MD10_20_08_p20.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.8 m40 >  logs/cluster_run_MD10_20_08_m40.log &                    
nohup python 004_6_clusters_simput.py MD10 20. 0.8 p40 >  logs/cluster_run_MD10_20_08_p40.log &                    

nohup python 004_6_clusters_simput.py MD10 2.  0.9 0 >    logs/cluster_run_MD10_2_09_0.log &                    
nohup python 004_6_clusters_simput.py MD10 2.  0.8 0 >    logs/cluster_run_MD10_2_08_0.log &                    
nohup python 004_6_clusters_simput.py MD10 2.  0.7 0 >    logs/cluster_run_MD10_2_07_0.log &                    

nohup python 004_6_clusters_simput.py MD40 2.  0.9 0 >    logs/cluster_run_MD40_2_09_0.log &                    
nohup python 004_6_clusters_simput.py MD40 2.  0.8 0 >    logs/cluster_run_MD40_2_08_0.log &                    
nohup python 004_6_clusters_simput.py MD40 2.  0.7 0 >    logs/cluster_run_MD40_2_07_0.log &                    

nohup python 004_6_clusters_simput.py MD40 20.  0.9 0 >   logs/cluster_run_MD40_20_09_0.log &                    
nohup python 004_6_clusters_simput.py MD40 20.  0.8 0 >   logs/cluster_run_MD40_20_08_0.log &                    
nohup python 004_6_clusters_simput.py MD40 20.  0.7 0 >   logs/cluster_run_MD40_20_07_0.log &                    

nohup python 004_6_clusters_simput.py MD04 20.  0.9 0 >   logs/cluster_run_MD04_2_09_0.log &                  
nohup python 004_6_clusters_simput.py MD04 20.  0.8 0 >   logs/cluster_run_MD04_2_08_0.log &                  
nohup python 004_6_clusters_simput.py MD04 20.  0.7 0 >   logs/cluster_run_MD04_2_07_0.log &                  

# computes cluster X-ray spectra
# depends on pyxspec
# script: 004_7_clusters_create_xpectra.py
# inputs: env, pyXspec
# outputs: ${env}/cat_CLU_SIMPUT/cluster_Xspectra/cluster_spectrum_10kT_????_100z_????.fits
# python 004_7_clusters_create_xpectra.py MD40
nohup nice -n 19 python 004_7_clusters_create_xpectra.py MD10 > run_CLU_7_MD10.log &
nohup nice -n 19 python 004_7_clusters_create_xpectra.py MD04 > run_CLU_7_MD04.log &
nohup nice -n 19 python 004_7_clusters_create_xpectra.py MD40 > run_CLU_7_MD40.log &
