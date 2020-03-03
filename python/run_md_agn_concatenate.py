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

stilts_cmd = 'stilts'

test_dir = os.path.join(os.environ[env], 'fits')
catalogue_dir = os.path.join(os.environ[env], 'cat_AGN_' + ftyp)

z_min = 0.0
if env == 'MD04':
    z_max = 0.45
else:
    z_min = 0.0
    z_max = 6.1

x_paint_git_dir = os.path.join(os.environ['GIT_XPAINT'])
lss_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python')

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


def run_all_agn(env, baseName):
	print('agn + eROSITA catalog', env, baseName, time.time() - t0)
	os.chdir(lss_git_dir)
	path_2_agn_file = os.path.join(test_dir, baseName + '_agn.fits')
	#path_2_eRO_catalog = os.path.join(test_dir, baseName + '_AGN.fit')
	dir_2_eRO_catalog = os.path.join(test_dir, 'cat_AGN_' + baseName)
	# AGN-related quantities
	if os.path.isfile(path_2_agn_file)==False:
		command_agn = "python3 003_0_agn.py " + env + ' ' + baseName
		print('=======================================')
		print(command_agn)
		print('=======================================')
		os.system(command_agn)
		# AGN fits catalog per pixel
		command_catalog = "python3 003_1_agn_catalogs.py " + env + ' ' + baseName
		print('=======================================')
		print(command_catalog)
		print('=======================================')
		os.system(command_catalog)

def run_agn_HPX_CAT(env, baseName):
	#print('agn + eROSITA catalog', env, baseName, time.time() - t0)
	os.chdir(lss_git_dir)
	path_2_agn_file = os.path.join(test_dir, baseName + '_agn.fits')
	#path_2_eRO_catalog = os.path.join(test_dir, baseName + '_AGN.fit')
	dir_2_eRO_catalog = os.path.join(test_dir, 'cat_AGN_' + baseName)
	# AGN-related quantities
	# AGN fits catalog per pixel
	command_catalog = "python3 003_1_agn_catalogs.py " + env + ' ' + baseName
	#print('=======================================')
	print('nohup ' + command_catalog + ' > 003_1_' + baseName+'.log & ')
	#print('=======================================')
	#os.system(command_catalog)


def concatenate_agn(env):
	if os.path.isdir(catalogue_dir) == False:
		os.system('mkdir -p ' + catalogue_dir)
	print('concatenate all AGN catalogs', env, time.time() - t0)
	for HEALPIX_32_id in n.arange(healpy.nside2npix(8)):
		path_2_eRO_catalogs = n.array(
			glob.glob(
				os.path.join(
					test_dir,
					'cat_AGN_' +
					ftyp +
					'_*',
					str(HEALPIX_32_id).zfill(6) +
					'.fit')))
		print(HEALPIX_32_id, len(path_2_eRO_catalogs), time.time() - t0)
		path_2_agn_summary_file = os.path.join(
			catalogue_dir, str(HEALPIX_32_id).zfill(6) + '.fit')
		print(path_2_agn_summary_file)
		os.chdir(lss_git_dir)
		# concatenates shells into a single fits catalog
		list_name = 'fit_list_' + env + "_" + ftyp + str(HEALPIX_32_id).zfill(6) + '_AGN.list'
		c1 = 'ls ' + test_dir + '/cat_AGN_' + ftyp + '_*/' + str(HEALPIX_32_id).zfill(6) + '.fit > ' + list_name
		print(c1)
		os.system(c1)
		c2 = stilts_cmd + """ tcat in=@""" + list_name + \
			""" ifmt=fits omode=out ofmt=fits out=""" + path_2_agn_summary_file
		print(c2)
		os.system(c2)
		os.system('rm ' + list_name)


#for bn in baseNames[::-1]:
	#print(env, bn)
	#run_all_agn(env, bn)
	##run_agn_HPX_CAT(env, bn)

concatenate_agn(env)

print('computes XLF')
command_plots = "python3 003_2_agn_compute_XLF_logNlogS_R.py " + env + ' ' + ftyp
print(command_plots)
os.system(command_plots)

print('plot logNlogS')
command_plots_2 = "python3 003_3_agn_plot_logNlogS.py " + env + ' ' + ftyp
print(command_plots_2)
os.system(command_plots_2)
