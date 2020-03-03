"""
commands for testing purpose :

rm $MD04/fits/*_galaxiesAroundClusters.fit
rm $MD10/fits/*_galaxiesAroundClusters.fit

python 004_1_galaxies_around_clusters.py MD10 all_0.95670 200c
python 004_1_galaxies_around_clusters.py MD10 all_0.89510 200c
python 004_1_galaxies_around_clusters.py MD04 all_0.95600 200c
python 004_1_galaxies_around_clusters.py MD04 all_0.88630 200c

python 004_2_cluster_galaxies.py MD10 all_0.95670 200c
python 004_3_cluster_red_galaxies.py MD10 all_0.95670 200c
python 004_4_red_sequence.py MD10 all_0.95670 200c

python 004_2_cluster_galaxies.py MD04 all_0.95600 200c
python 004_3_cluster_red_galaxies.py MD04 all_0.95600 200c
python 004_4_red_sequence.py MD04 all_0.95600 200c

python 004_2_cluster_galaxies.py MD10 all_0.89510 200c
python 004_3_cluster_red_galaxies.py MD10 all_0.89510 200c
python 004_4_red_sequence.py MD10 all_0.89510 200c

python 004_2_cluster_galaxies.py MD04 all_0.88630 200c
python 004_3_cluster_red_galaxies.py MD04 all_0.88630 200c
python 004_4_red_sequence.py MD04 all_0.88630 200c

python 004_5_plot_clusters.py MD04 200c
python 004_5_plot_clusters.py MD10 200c

python 004_5_plot_clusters_colors.py MD04 
python 004_5_plot_clusters_colors.py MD10

python 004_5_plot_clusters_richness.py MD10 200c
python 004_5_plot_clusters_richness.py MD04 200c
"""
import os, sys

env =  sys.argv[1] # 'MD04'
baseName = sys.argv[2]  # "all_0.62840"
delta_crit = sys.argv[3] # '200c'
print(env, baseName, delta_crit)

os.system("python3 004_2_cluster_galaxies.py "+env+" "+baseName+" "+delta_crit)
os.system("python3 004_3_cluster_red_galaxies.py "+env+" "+baseName+" "+delta_crit)
os.system("python3 004_4_red_sequence.py "+env+" "+baseName+" "+delta_crit)





