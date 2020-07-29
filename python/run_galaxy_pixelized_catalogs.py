import glob
import sys
from astropy_healpix import healpy
import os
import numpy as n
import time
t0 = time.time()

env = sys.argv[1] 
ftyp = 'all' # or 'sat'
test_dir = os.path.join(os.environ[env], 'fits')
catalogue_dir = os.path.join(os.environ[env], 'cat_GALAXY_' + ftyp)

z_min = 0.0
z_max = 6.1

lss_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python')

baseNames_all = sorted(n.array([os.path.basename(el)[:-12] for el in n.array(glob.glob(os.path.join(test_dir, ftyp +'_?.?????_galaxy.fits')))]))

z_array = n.array([1. / float(eel.split('_')[1]) - 1 for eel in baseNames_all])
z_sel = (z_array < z_max) & (z_array >= z_min)
baseNames = n.array(baseNames_all)[z_sel]
print(baseNames)


def run_pixelize_galaxies(env, baseName):
	command_catalog = "python3 002_1_galaxy_catalogs.py " + env + ' ' + baseName
	#print('=======================================')
	#print(command_catalog)
	print('nohup '+ command_catalog + '> logs/002_1_galCAT_'+env + "_" + baseName + '.log &')
	dir_2_eRO_catalog = os.path.join(test_dir, 'cat_GALAXY_' + baseName)
	#if os.path.isdir(dir_2_eRO_catalog) == False:
		#os.system(command_catalog)
	#else:
		#print('already done')

def concatenate_galaxies(env):
	if os.path.isdir(catalogue_dir) == False:
		os.system('mkdir -p ' + catalogue_dir)
	for HEALPIX_id in n.arange(healpy.nside2npix(8)):
		path_2_galaxy_summary_file = os.path.join(catalogue_dir, str(HEALPIX_id).zfill(6) + '.fit')
		#print(path_2_galaxy_summary_file)
		os.chdir(lss_git_dir)
		# concatenates shells into a single fits catalog
		list_name = 'logs/fit_list_' + env + "_" + ftyp + str(HEALPIX_id).zfill(6) + '_GALAXY.list'
		c1 = 'ls ' + test_dir + '/cat_GALAXY_' + ftyp + '_*/' + str(HEALPIX_id).zfill(6) + '.fit > ' + list_name
		print(c1)
		#os.system(c1)
		c2 = 'python concat_file_list.py '+list_name+' '+path_2_galaxy_summary_file
		command_2 = 'nohup '+c2 + '> logs/002_1_concat_'+env + "_" + ftyp + str(HEALPIX_id).zfill(6)+'.log &'
		#print(command_2)
		#os.system(c2)
		#os.system('rm ' + list_name)

def concatenate_galaxies_c2(env):
	if os.path.isdir(catalogue_dir) == False:
		os.system('mkdir -p ' + catalogue_dir)
	for HEALPIX_id in n.arange(healpy.nside2npix(8)):
		path_2_galaxy_summary_file = os.path.join(catalogue_dir, str(HEALPIX_id).zfill(6) + '.fit')
		#print(path_2_galaxy_summary_file)
		os.chdir(lss_git_dir)
		# concatenates shells into a single fits catalog
		list_name = 'logs/fit_list_' + env + "_" + ftyp + str(HEALPIX_id).zfill(6) + '_GALAXY.list'
		c1 = 'ls ' + test_dir + '/cat_GALAXY_' + ftyp + '_*/' + str(HEALPIX_id).zfill(6) + '.fit > ' + list_name
		#print(c1)
		#os.system(c1)
		c2 = 'python concat_file_list.py '+list_name+' '+path_2_galaxy_summary_file
		command_2 = 'nohup '+c2 + '> logs/002_1_concat_'+env + "_" + ftyp + str(HEALPIX_id).zfill(6)+'.log &'
		command_2 = c2 # + '> logs/002_1_concat_'+env + "_" + ftyp + str(HEALPIX_id).zfill(6)+'.log &'
		print(command_2)
		#os.system(c2)
		#os.system('rm ' + list_name)

print('#============================================================')
print('#============================================================')
print('# pixelization ')
print('#============================================================')
print('#============================================================')

for bn in baseNames[::-1]:
	#print(env, bn)
	run_pixelize_galaxies(env, bn)

print('#============================================================')
print('#============================================================')
print('# concatenate ')
print('#============================================================')
print('#============================================================')

concatenate_galaxies(env)
concatenate_galaxies_c2(env)
