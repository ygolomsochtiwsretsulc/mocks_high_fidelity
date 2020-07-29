import glob
import sys
import os
import numpy as n
import time
t0 = time.time()

env =  sys.argv[1]  # 'MD10'
ftyp = 'all' # or 'sat'

test_dir = os.path.join(os.environ[env], 'fits')

z_min = 0.0
z_max = 6.1

lss_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python')
lss_fig_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'figures')

baseNames_all = sorted(
    n.array(
        [
            os.path.basename(el)[
                :-
                5] for el in n.array(
                    glob.glob(
                        os.path.join(
                            test_dir,
                            ftyp +
                            '_?.?????.fits')))]))

z_array = n.array([1. / float(eel.split('_')[1]) - 1 for eel in baseNames_all])
z_sel = (z_array < z_max) & (z_array >= z_min)
baseNames = n.array(baseNames_all)[z_sel]
print(baseNames)


def run_all_coord(env, baseName):
	print('coordinates', env, baseName, time.time() - t0)
	os.chdir(lss_git_dir)
	# Adds sky coordinates + redshift distances
	command = "python3 001_coordinates.py " + env + ' ' + baseName
	print(command)
	os.system(command)


def run_all_gal(env, baseName):
	path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
	path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
	os.chdir(lss_git_dir)
	if os.path.isfile(path_2_coordinate_file)==False:
		#print('coordinates + galaxies', env, baseName, time.time() - t0)
		# Adds sky coordinates + redshift distances
		command = "python3 001_coordinates.py " + env + ' ' + baseName
		print("nohup "+command+" > logs/"+env+'_'+baseName+'001_0_coord.log'+" &")
		os.system(command)
	if os.path.isfile(path_2_galaxy_file)==False:
		# Adds galaxy properties
		command = "python3 002_0_galaxy.py " + env + ' ' + baseName
		print("nohup "+command+" > logs/"+env+'_'+baseName+'002_0_galaxy.log'+" &")
		os.system(command)

def run_both(env, baseName):
	os.chdir(lss_git_dir)
	print('coordinates + galaxies', env, baseName, time.time() - t0)
	# Adds sky coordinates + redshift distances
	command = "python3 001_coordinates.py " + env + ' ' + baseName
	print(command)
	os.system(command)		# Adds galaxy properties
	command = "python3 002_0_galaxy.py " + env + ' ' + baseName
	print(command)
	os.system(command)

def run_both_galaxy_v0(env, baseName):
	os.chdir(lss_git_dir)
	print('coordinates + galaxies', env, baseName, time.time() - t0)
	# Adds sky coordinates + redshift distances
	command = "python3 001_coordinates.py " + env + ' ' + baseName
	print(command)
	os.system(command)		# Adds galaxy properties
	command = "python3 002_0_galaxy_v0.py " + env + ' ' + baseName
	print(command)
	os.system(command)
	
def run_gal_only(env, baseName):
	path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
	os.chdir(lss_git_dir)
	print('galaxies', env, baseName, time.time() - t0)
	command = "python3 002_0_galaxy.py " + env + ' ' + baseName
	print(command)
	os.system(command)

def plot_stellar_mass(env, baseName):
	os.chdir(lss_fig_git_dir)
	print('galaxies', env, baseName, time.time() - t0)
	command = "python3 plot_SMHMR.py " + env + ' ' + baseName
	print(command)
	os.system(command)

for bn in baseNames[::-1]:
	#print(env, bn)
	run_both_galaxy_v0(env, bn)
	#run_all_gal(env, bn)
	#run_gal_only(env, bn)
	#plot_stellar_mass(env, bn)
