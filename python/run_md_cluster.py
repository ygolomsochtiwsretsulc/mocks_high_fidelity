import glob
import sys
import os
import numpy as n
import time
t0 = time.time()

env = sys.argv[1]  # 'MD10'
with_image = sys.argv[2]  # 'yes'

test_dir = os.path.join(os.environ[env], 'fits')

baseNames_all = n.array([os.path.basename(el)[:-5] for el in n.array(glob.glob(os.path.join(test_dir, 'all_?.?????.fits')))])
baseNames_all.sort()
print(baseNames_all)
if env == 'MD04':
	z_max = 0.45
else:
	z_max = 1.2

z_array = n.array([1. / float(eel.split('_')[1]) - 1 for eel in baseNames_all])
z_sel = (z_array < z_max)
baseNames = baseNames_all[z_sel]
baseNames.sort()
print(baseNames)


def run_all_clusters(env, baseName):
	print('cluster', env, baseName, time.time() - t0)
	# pathes to cluster output files :
	# creates input for X-ray painting code
	path_2_CLU_file = os.path.join(test_dir, baseName + '_CLU.fits')
	print('mCLU file', path_2_CLU_file, time.time() - t0)
	# parent cluster catalog M500c>5e13
	#if os.path.isfile(path_2_CLU_file)==False:
	os.system("python3 004_0_cluster.py " + env + ' ' + baseName + ' ' + with_image)


for bn in baseNames[::-1]:
	print(env, bn)
	run_all_clusters(env, bn)

#os.system("python3 004_1_cluster_Merge.py " + env )
#os.system("python3 004_2_cluster_galaxies.py " + env )
#os.system("python3 004_3_cluster_red_galaxies.py " + env)
#os.system("python3 004_4_red_sequence.py " + env)
