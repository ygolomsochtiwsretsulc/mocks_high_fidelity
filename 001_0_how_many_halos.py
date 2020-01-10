"""
What it does
------------

counts the number of halos 


"""
# import python packages

import sys, os, glob
import astropy.io.fits as fits
from astropy.table import Table, Column
import numpy as n
import time
print('count halos')

# import all pathes
envs = n.array(['MD10', 'MD04', 'MD40', 'UNIT_fA1i_DIR', 'UNIT_fA1_DIR', 'UNIT_fA2_DIR'])
f = open('halo_count.txt', 'w')

for env in envs:
	t0 = time.time()
	f.write(env)
	test_dir = os.path.join(os.environ[env], 'fits')
	path_2_light_cones = glob.glob(os.path.join(test_dir, 'all_?.?????.fits'))
	path_2_light_cones.sort()
	NN = []
	for path_2_light_cone in path_2_light_cones[::-1]:
		f1 = fits.open(path_2_light_cone)
		PIDS = f1[1].data['pid']
		distinct = (PIDS==-1)
		N_all = len(PIDS)
		N_d = len(PIDS[distinct])
		f.write('   == '+ os.path.basename(path_2_light_cone)[4:-5]+' '+str( N_d)+' '+str( N_all-N_d))
		f1.close()
		NN.append([N_all, N_d])
	NN = n.transpose(NN)
	f.write(env+' '+str( NN.sum(axis=1)[1])+ ' & '+str( NN.sum(axis=1)[0]-NN.sum(axis=1)[1])+'\n' )
	#print('+++++++++++++++++++++++++++++++')
	#print('+++++++++++++++++++++++++++++++')
	#print('+++++++++++++++++++++++++++++++')
f.close()