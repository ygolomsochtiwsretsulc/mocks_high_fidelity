"""
Creates a fits catalog containing the 4FS input columns.

Create the ID array
remove unwanted column
check rulesets and templates names are correct

"""
import glob
import sys
from astropy_healpix import healpy
import os
from scipy.stats import scoreatpercentile
import pandas as pd  # external package
from scipy.special import erf
from astropy.coordinates import SkyCoord
import astropy.constants as cc
import astropy.io.fits as fits
from astropy.table import Table, Column
import astropy.units as u
import numpy as n
import extinction
print('CREATES 4FS FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')

# astropy
#from lib_model_agn import *
#option = sys.argv[1]
#option = 'SNR3'
#option = 'eRASS8'
#option = 'eRASS3'

env = sys.argv[1]

root_dir = os.path.join(os.environ[env])

dir_2_eRO_all = os.path.join(root_dir, "cat_AGN-MAG_all")

dir_2_4MOST = os.path.join(root_dir, "cat_AGN_4MOST")

if os.path.isdir(dir_2_4MOST) == False:
    os.system('mkdir -p ' + dir_2_4MOST)

path_2_flux_limits = os.path.join(os.environ['GIT_AGN_MOCK'], "data", "erosita", "flux_limits.fits")
flux_lim_data = fits.open(path_2_flux_limits) 
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3
FX_LIM = flux_lim_data[1].data['flux_limit_SNR3']
#FX_LIM = n.log10(2e-17 )


area_per_cat = healpy.nside2pixarea(8, degrees=True)
N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels):
	#HEALPIX_id = 350
	print(HEALPIX_id)
	path_2_eRO_all_catalog = os.path.join(dir_2_eRO_all, str(HEALPIX_id).zfill(6) + '.fit')
	path_2_AGN_WIDE_catalog = os.path.join( dir_2_4MOST, 'AGN_WIDE_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_AGN_DEEP_catalog = os.path.join( dir_2_4MOST, 'AGN_DEEP_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_AGN_IR_catalog = os.path.join( dir_2_4MOST, 'AGN_IR_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_QSO_catalog = os.path.join( dir_2_4MOST, 'QSO_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_LyA_catalog = os.path.join( dir_2_4MOST, 'LyA_' + str(HEALPIX_id).zfill(6) + '.fit')

	bitlist = {'AGN_WIDE':0, 'AGN_DEEP':1, 'AGN_IR':2, 'QSO':3, 'LyA':4 }
	priority_values = {'AGN_WIDE':100, 'AGN_DEEP':100, 'AGN_IR':99, 'QSO':99, 'LyA':100 }

	# reassigns templates correctly
	z_all = n.hstack(( 0., n.arange(0.3, 3., 0.2), 3.5, 4.5, 6. ))
	zmins = z_all[:-1]
	zmaxs = z_all[1:]

	hd_all = fits.open(path_2_eRO_all_catalog)
	N_agn_all = len(hd_all[1].data['ra'])

	pix_ids = healpy.ang2pix(512, n.pi/2. - hd_all[1].data['g_lat'] *n.pi /180., hd_all[1].data['g_lon']*n.pi /180., nest=True)
	FX_LIM_value_cen = 10**(FX_LIM[pix_ids] - 0.1 )
	# =============================
	# =============================
	# eROSITA SAMPLE
	# apply flux limit erosita
	# =============================
	# =============================
	# X-ray
	detected_all = (hd_all[1].data['AGN_FX_soft'] > FX_LIM_value_cen) & ( hd_all[1].data['g_lon'] > 180)
	print(len((detected_all).nonzero()[0])*768 )
	# optical cut
	opt_all = ( hd_all[1].data['AGN_SDSS_r_magnitude'] < 23.5 )
	# overall erosita selections
	ero_all = ( detected_all ) & ( opt_all )
	# optical t1
	opt_t11_all = (opt_all) & (hd_all[1].data['AGN_type']==11)
	# type 2 AGN
	T2_all = (hd_all[1].data['AGN_type'] == 22) | (hd_all[1].data['AGN_type'] == 21)
	# wise selection
	wise_def_all = (hd_all[1].data['WISE-W1'] < 26) & (hd_all[1].data['WISE-W1_err'] < 1) & (hd_all[1].data['WISE-W1'] < 19 ) # 18.85)
	# IR sample, not detected in X-ray i.e. additional targets
	agn_ir_all =   (T2_all) & (opt_all) & (wise_def_all) # ( detected_all == False) &
	# area selection
	area_erosita_all = (abs(hd_all[1].data['g_lat']) > 10) # #& (hd_all[1].data['dec'] < 10)
	# area opt-IR
	area_optIR_all = (abs(hd_all[1].data['g_lat']) > 10) #& (hd_all[1].data['dec'] < 10)
	# selection booleans
	is_erosita_all = ( (ero_all) & (area_erosita_all) )
	is_optIR_all   = ( (agn_ir_all) & (area_optIR_all) ) 
	is_optT11_all  = ( (opt_t11_all) & (area_optIR_all) )
	# complete selection for all 4most sub surveys
	sel_all = ( is_erosita_all ) | ( is_optIR_all ) | ( is_optT11_all )
	# number per square degrees
		# =============================
	# Define sub surveys
	# =============================
	# initialize target bits

	survbit_all = n.zeros(len(hd_all[1].data['dec']),dtype='int')
	# subsurvey selections 
	# S6, AGN_DEEP
	AGN_DEEP_all = ( is_erosita_all ) & ( abs(hd_all[1].data['ecl_lat']) > 80 )
	# S6, AGN_WIDE
	AGN_WIDE_all = ( is_erosita_all ) & (AGN_DEEP_all==False) # ( abs(hd_all[1].data['ecl_lat']) < 80 )
	# S6, AGN_IR 
	AGN_IR_all = ( is_optIR_all ) 
	# S8, QSO
	QSO_all = ( is_optT11_all ) & ( hd_all[1].data['redshift_R'] < 2.2 )
	# S8, Lya
	LyA_all = ( is_optT11_all ) & ( hd_all[1].data['redshift_R'] >= 2.2 ) 
	# assign target bits
	survbit_all[AGN_WIDE_all] += 2**bitlist['AGN_WIDE']
	#
	survbit_all[AGN_DEEP_all] += 2**bitlist['AGN_DEEP']
	#
	survbit_all[AGN_IR_all] += 2**bitlist['AGN_IR']
	#
	survbit_all[QSO_all] += 2**bitlist['QSO']
	#
	survbit_all[LyA_all] += 2**bitlist['LyA']
	#

	target_bit = survbit_all[sel_all]


	# Now create a 4MOST file per sub survey
	# merge cen + sat
	def get_table(subsurvey):
		t_survey = Table()
		bit_value = bitlist[subsurvey]
		sel = ( target_bit & 2**bit_value != 0)

		def create_column(col_name):
			return hd_all[1].data[col_name][sel_all][sel]

		ra_array = create_column('ra')
		dec_array = create_column('dec')
		N_obj = len(ra_array)
		if N_obj>0:
			N1 = n.arange(N_obj)
			id_list = HEALPIX_id*1e8 + N1
			NAME = n.array([ str(int(el)).zfill(11) for el in id_list ])
			t_survey.add_column(Column(name='NAME', data=NAME, unit=''))
			t_survey.add_column(Column(name='RA', data=ra_array, unit='deg'))
			t_survey.add_column(Column(name='DEC', data=dec_array, unit='deg'))
			PMRA = n.zeros(N_obj)
			t_survey.add_column(Column(name='PMRA', data=PMRA, unit='mas/yr'))
			PMDEC = n.zeros(N_obj)
			t_survey.add_column(Column(name='PMDEC', data=PMDEC, unit='mas/yr'))
			EPOCH = n.ones(N_obj)*2000.
			t_survey.add_column(Column(name='EPOCH', data=EPOCH, unit='yr'))
			# 'RESOLUTION':n.int16, 1I
			RESOLUTION = n.ones(N_obj).astype('int')
			t_survey.add_column(Column(name='RESOLUTION', data=RESOLUTION, unit=''))

			SUBSURVEY = n.ones(N_obj).astype('str')
			SUBSURVEY[:] = subsurvey
			t_survey.add_column(Column(name='SUBSURVEY', data=SUBSURVEY, unit=''))
			# 'PRIORITY':n.int16, 1I
			PRIORITY = n.zeros(N_obj).astype('int') + priority_values[subsurvey]
			t_survey.add_column(Column(name='PRIORITY', data=PRIORITY, unit=''))

			galactic_ebv_array = create_column('galactic_ebv')
			t_survey.add_column(Column(name='REDDENING',data=galactic_ebv_array, unit='mag'))

			# REDDENING for templates
			ebv_1000 = (t_survey['REDDENING']*1000).astype('int')
			#print('EBV', n.min(ebv_1000), n.max(ebv_1000))
			ebv_1_0 = ( ebv_1000 > 1000 ) 
			ebv_0_5 = ( ebv_1000 > 500 ) & ( ebv_1000 <= 1000 ) 
			ebv_0_4 = ( ebv_1000 > 400 ) & ( ebv_1000 <= 500 ) 
			ebv_0_3 = ( ebv_1000 > 300 ) & ( ebv_1000 <= 400 ) 
			ebv_0_2 = ( ebv_1000 > 200 ) & ( ebv_1000 <= 300 ) 
			ebv_0_1 = ( ebv_1000 > 100 ) & ( ebv_1000 <= 200 ) 
			ebv_0_0 = ( ebv_1000 <= 100 ) 
			z_name = lambda z0, z1 : "_zmin_"+str(int(10*z0)).zfill(2)+"_zmax_"+str(int(10*z1)).zfill(2)
			# templates
			template_names = n.zeros(N_obj).astype('U100')
			#
			ruleset_array = n.zeros(N_obj).astype('str')
			# 
			if subsurvey == "AGN_WIDE" or subsurvey == "AGN_DEEP" or subsurvey == "AGN_IR":
				ruleset_array[:] = "AGN_ALL_3PC"
			if subsurvey == "QSO" or subsurvey == "LyA":
				ruleset_array[:] = "COSMO_AGN"

			AGN_type_array = create_column('AGN_type')
			AGN_random_number_array = create_column('AGN_random_number')
			z_array = create_column('redshift_R')

			QSO = (AGN_type_array == 11) | (AGN_type_array == 12)
			T2 = (AGN_type_array == 22) | (AGN_type_array == 21)
			ELL = (T2) & (AGN_random_number_array < 0.2)

			for z0, z1 in zip(zmins, zmaxs):
				zsel = (z_array >= z0) & (z_array < z1)
				template_names[(zsel)] = "4most_" + 'qso_BL' +  z_name(z0, z1) + '_EBV_0_01.fits'
				template_names[(QSO) & (zsel) & (ebv_0_0)] = "4most_" +  'qso_BL' + z_name(z0, z1) + '_EBV_0_01.fits'
				template_names[(QSO) & (zsel) & (ebv_0_1)] = "4most_" + 'qso_BL' + z_name(z0, z1) + '_EBV_0_1.fits'
				template_names[(QSO) & (zsel) & (ebv_0_2)] = "4most_" + 'qso_BL' + z_name(z0, z1) + '_EBV_0_2.fits'
				template_names[(QSO) & (zsel) & (ebv_0_3)] = "4most_" + 'qso_BL' + z_name(z0, z1) + '_EBV_0_3.fits'
				template_names[(QSO) & (zsel) & (ebv_0_4)] = "4most_" + 'qso_BL' + z_name(z0, z1) + '_EBV_0_4.fits'
				template_names[(QSO) & (zsel) & (ebv_0_5)] = "4most_" + 'qso_BL' + z_name(z0, z1) + '_EBV_0_5.fits'
				template_names[(QSO) & (zsel) & (ebv_1_0)] = "4most_" + 'qso_BL' + z_name(z0, z1) + '_EBV_1_0.fits'
				if z1 < 2.2:
					template_names[(T2) & (zsel) & (ebv_0_0)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
					template_names[(T2) & (zsel) & (ebv_0_1)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
					template_names[(T2) & (zsel) & (ebv_0_2)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
					template_names[(T2) & (zsel) & (ebv_0_3)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
					template_names[(T2) & (zsel) & (ebv_0_4)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
					template_names[(T2) & (zsel) & (ebv_0_5)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
					template_names[(T2) & (zsel) & (ebv_1_0)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'

					template_names[(ELL) & (zsel) & (ebv_0_0)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_0_01.fits'
					template_names[(ELL) & (zsel) & (ebv_0_1)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_0_1.fits'
					template_names[(ELL) & (zsel) & (ebv_0_2)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_0_2.fits'
					template_names[(ELL) & (zsel) & (ebv_0_3)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_0_3.fits'
					template_names[(ELL) & (zsel) & (ebv_0_4)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_0_4.fits'
					template_names[(ELL) & (zsel) & (ebv_0_5)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_0_5.fits'
					template_names[(ELL) & (zsel) & (ebv_1_0)] = "4most_" + 'LRG' + z_name(z0, z1) + '_EBV_1_0.fits'
				if z1 >= 2.2 and z1 < 6.:
					template_names[(T2) & (zsel) & (ebv_0_0)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
					template_names[(T2) & (zsel) & (ebv_0_1)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
					template_names[(T2) & (zsel) & (ebv_0_2)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
					template_names[(T2) & (zsel) & (ebv_0_3)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
					template_names[(T2) & (zsel) & (ebv_0_4)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
					template_names[(T2) & (zsel) & (ebv_0_5)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
					template_names[(T2) & (zsel) & (ebv_1_0)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'
					template_names[(ELL) & (zsel) & (ebv_0_0)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
					template_names[(ELL) & (zsel) & (ebv_0_1)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
					template_names[(ELL) & (zsel) & (ebv_0_2)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
					template_names[(ELL) & (zsel) & (ebv_0_3)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
					template_names[(ELL) & (zsel) & (ebv_0_4)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
					template_names[(ELL) & (zsel) & (ebv_0_5)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
					template_names[(ELL) & (zsel) & (ebv_1_0)] = "4most_" + 'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'

			# 'TEMPLATE':str, max 256 char
			t_survey.add_column(Column(name='TEMPLATE', data=template_names, unit=''))
			# 'RULESET':str, max 256 char
			t_survey.add_column(Column(name='RULESET', data=ruleset_array, unit=''))
			# 'REDSHIFT_ESTIMATE':n.float32, 1E
			# 'REDSHIFT_ERROR':n.float32, 1E
			t_survey.add_column(Column(name='REDSHIFT_ESTIMATE', data=z_array, unit=''))
			t_survey.add_column(Column(name='REDSHIFT_ERROR', data=n.ones(N_obj), unit=''))
			# 'EXTENT_FLAG': 1I
			# =1
			# 'EXTENT_PARAMETER': 1E
			# =0
			# 'EXTENT_INDEX': 1E
			# =0
			t_survey.add_column(Column(name='EXTENT_FLAG'     , data=n.zeros(N_obj).astype('int') , unit=''))
			t_survey.add_column(Column(name='EXTENT_PARAMETER', data=n.zeros(N_obj), unit=''))
			t_survey.add_column(Column(name='EXTENT_INDEX'    , data=n.zeros(N_obj), unit=''))
			# 'MAG':n.float32,
			# 'MAG_ERR':n.float32
			# 'MAG_TYPE': str max 256 char
			HSC_RMAG_array = create_column('HSC-r')
			HSC_RMAGERR_array = create_column('HSC-r_err')
			#
			r_v=3.1
			a_v = galactic_ebv_array * r_v
			delta_mag = n.hstack(( n.array([ extinction.fitzpatrick99(n.array([6500.]), el, r_v=3.1, unit='aa') for el in a_v ]) ))
			#rv = av/ebv
			#av = rv x ebv
			extincted_mag = HSC_RMAG_array + delta_mag
			t_survey.add_column(Column(name='MAG', data=extincted_mag, unit='mag'))
			t_survey.add_column(Column(name='MAG_ERR', data=HSC_RMAGERR_array, unit='mag'))
			MAG_TYPE = n.ones(N_obj).astype('str')
			MAG_TYPE[:] = 'DECam_r_AB'
			t_survey.add_column(Column(name='MAG_TYPE', data=MAG_TYPE, unit=''))
			# 'REDDENING':n.float32, 1E
			# 'DATE_EARLIEST':n.float64, JulianDate decimal days # 01-Nov-2022
			# 'DATE_LATEST':n.float64, JulianDate decimal days # 02-Feb-2033
			t_survey.add_column(Column(name='DATE_EARLIEST',data=22305 * n.ones(N_obj), unit='d'))
			t_survey.add_column(Column(name='DATE_LATEST'  ,data=33033 * n.ones(N_obj), unit='d'))

			t_survey.add_column(Column(name='target_bit'  ,data=target_bit[sel], unit=''))

			t_survey.add_column(Column(name='dL_cm'  ,data=create_column('dL_cm'), unit='cm'))
			t_survey.add_column(Column(name='galactic_NH'  ,data=create_column('galactic_NH'), unit=''))
			t_survey.add_column(Column(name='galaxy_stellar_mass'  ,data=create_column('galaxy_stellar_mass'), unit=''))
			t_survey.add_column(Column(name='HALO_Mvir'  ,data=create_column('HALO_Mvir'), unit=''))
			t_survey.add_column(Column(name='AGN_LX_soft'  ,data=create_column('AGN_LX_soft'), unit=''))
			t_survey.add_column(Column(name='AGN_FX_soft'  ,data=create_column('AGN_FX_soft'), unit=''))
			t_survey.add_column(Column(name='AGN_LX_hard'  ,data=create_column('AGN_LX_hard'), unit=''))
			t_survey.add_column(Column(name='AGN_FX_hard'  ,data=create_column('AGN_FX_hard'), unit=''))
			t_survey.add_column(Column(name='AGN_nH'  ,data=create_column('AGN_nH'), unit=''))
			t_survey.add_column(Column(name='WISE-W1'  ,data=create_column('WISE-W1'), unit=''))
			t_survey.add_column(Column(name='AGN_SDSS_r_magnitude'  ,data=create_column('AGN_SDSS_r_magnitude'), unit=''))
			t_survey.add_column(Column(name='g_lat'  ,data=create_column('g_lat'), unit='deg'))
			t_survey.add_column(Column(name='g_lon'  ,data=create_column('g_lon'), unit='deg'))
			t_survey.add_column(Column(name='ecl_lat'  ,data=create_column('ecl_lat'), unit='deg'))
			#t_survey.add_column(Column(name=''  ,data=create_column(''), unit=''))
			return t_survey
		else:
			return 0.

	subsurvey = 'AGN_WIDE'
	t_out = get_table(subsurvey)
	if t_out != 0:
		t_out.write (path_2_AGN_WIDE_catalog  , overwrite=True)

	subsurvey = 'AGN_DEEP'
	t_out = get_table(subsurvey)
	if t_out != 0:
		t_out.write (path_2_AGN_DEEP_catalog  , overwrite=True)

	subsurvey = 'AGN_IR'
	t_out = get_table(subsurvey)
	if t_out != 0:
		t_out.write (path_2_AGN_IR_catalog  , overwrite=True)

	subsurvey = 'QSO'
	t_out = get_table(subsurvey)
	if t_out != 0:
		t_out.write (path_2_QSO_catalog  , overwrite=True)

	subsurvey = 'LyA'
	t_out = get_table(subsurvey)
	if t_out != 0:
		t_out.write (path_2_LyA_catalog  , overwrite=True)

