import glob
import sys
import os
import numpy as n
import time
t0 = time.time()

env = sys.argv[1]  
with_image = sys.argv[2] 
pixel_size_image = sys.argv[3]
b_HS = sys.argv[4]
cov_mat_option = sys.argv[5] 
z_max = float( sys.argv[6] )

baseNames_all = n.array([os.path.basename(el)[:-5] for el in n.array(glob.glob(os.path.join(os.environ[env], 'fits', 'all_?.?????.fits')))])
baseNames_all.sort()

z_array = n.array([1. / float(eel.split('_')[1]) - 1 for eel in baseNames_all])
z_sel = (z_array < z_max)
baseNames = baseNames_all[z_sel]
baseNames.sort()
print(baseNames)

def run_all_clusters(env, baseName):
	command = "python 004_0_cluster.py " + env + ' ' + baseName + ' ' + with_image + ' ' + pixel_size_image + ' ' + b_HS + ' ' + cov_mat_option
	print(command)
	os.system(command)

for bn in baseNames[::-1]:
	run_all_clusters(env, bn)
