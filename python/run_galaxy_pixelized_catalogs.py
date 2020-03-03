"""
export GIT_AGN_MOCK='/home/comparat/software/lss_mock_dev/'
export GIT_XPAINT='/home/comparat/software/XSB_Painting'


"""
import glob
import sys
from astropy_healpix import healpy
import os
import numpy as n
import time
t0 = time.time()

env = sys.argv[1]  # 'MD10'
ftyp = sys.argv[2]  # 'all' or 'sat'
#laptop = sys.argv[3]  # 'True'

#if laptop == "True":
    #stilts_cmd = 'java -jar /home/comparat/software/stilts.jar'
#else:
stilts_cmd = 'stilts'

test_dir = os.path.join(os.environ[env], 'fits')
catalogue_dir = os.path.join(os.environ[env], 'cat_GALAXY_' + ftyp)

z_min = 0.0
if env == 'MD04':
	z_max = 0.45
else:
	z_min = 0.
	z_max = 6.1

lss_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python')

baseNames_all = sorted(n.array([os.path.basename(el)[:-5] for el in n.array(glob.glob(os.path.join(test_dir, ftyp +'_?.?????.fits')))]))

z_array = n.array([1. / float(eel.split('_')[1]) - 1 for eel in baseNames_all])
z_sel = (z_array < z_max) & (z_array >= z_min)
baseNames = n.array(baseNames_all)[z_sel]
print(baseNames)


def run_pixelize_galaxies(env, baseName):
	command_catalog = "python3 002_1_galaxy_catalogs.py " + env + ' ' + baseName
	print('=======================================')
	print(command_catalog)
	dir_2_eRO_catalog = os.path.join(test_dir, 'cat_GALAXY_' + baseName)
	if os.path.isdir(dir_2_eRO_catalog) == False:
		os.system(command_catalog)
	else:
		print('already done')

def concatenate_galaxies(env):
	if os.path.isdir(catalogue_dir) == False:
		os.system('mkdir -p ' + catalogue_dir)
	print('concatenate all GALAXY catalogs', env, time.time() - t0)
	for HEALPIX_id in n.arange(healpy.nside2npix(8)):
		path_2_galaxy_summary_file = os.path.join(catalogue_dir, str(HEALPIX_id).zfill(6) + '.fit')
		print(path_2_galaxy_summary_file)
		os.chdir(lss_git_dir)
		# concatenates shells into a single fits catalog
		list_name = 'fit_list_' + env + "_" + ftyp + str(HEALPIX_id).zfill(6) + '_GALAXY.list'
		c1 = 'ls ' + test_dir + '/cat_GALAXY_' + ftyp + '_*/' + str(HEALPIX_id).zfill(6) + '.fit > ' + list_name
		print(c1)
		os.system(c1)
		c2 = stilts_cmd + """ tcat in=@""" + list_name + """ ifmt=fits omode=out ofmt=fits out=""" + path_2_galaxy_summary_file
		print(c2)
		os.system(c2)
		os.system('rm ' + list_name)


for bn in baseNames[::-1]:
	print(env, bn)
	run_pixelize_galaxies(env, bn)

concatenate_galaxies(env)
