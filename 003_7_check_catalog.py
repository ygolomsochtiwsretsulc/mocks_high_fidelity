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
print('VERIFICATION OF 4FS FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')

env = 'MD10'
root_dir = os.path.join(os.environ[env])

nl = lambda selection : len(selection.nonzero()[0])

dir_2_AGN_DEEP = os.path.join(root_dir, "AGN_DEEP_4MOST.fits")
dir_2_AGN_WIDE = os.path.join(root_dir, "AGN_WIDE_4MOST.fits")
dir_2_AGN_IR   = os.path.join(root_dir, "AGN_IR_4MOST.fits")
dir_2_QSO      = os.path.join(root_dir, "QSO_4MOST.fits")
dir_2_LyA      = os.path.join(root_dir, "LyA_4MOST.fits")

t_AGN_DEEP = Table.read(dir_2_AGN_DEEP )
t_AGN_WIDE = Table.read(dir_2_AGN_WIDE )
t_AGN_IR   = Table.read(dir_2_AGN_IR   )
t_QSO      = Table.read(dir_2_QSO      )
t_LyA      = Table.read(dir_2_LyA      )

bitlist = {'AGN_WIDE':0, 'AGN_DEEP':1, 'AGN_IR':2, 'QSO':3, 'LyA':4 }

def compute_stat(t_survey = t_AGN_DEEP):
	s0 = (t_survey['target_bit'] & 2**bitlist['AGN_WIDE'] != 0 )
	s1 = (t_survey['target_bit'] & 2**bitlist['AGN_DEEP'] != 0 )
	s2 = (t_survey['target_bit'] & 2**bitlist['AGN_IR']   != 0 )
	s3 = (t_survey['target_bit'] & 2**bitlist['QSO']      != 0 )
	s4 = (t_survey['target_bit'] & 2**bitlist['LyA']      != 0 )
	#print(nl(s0), nl(s1), nl(s2), nl(s3), nl(s4))
	return str(nl(s0)), str(nl(s1)), str(nl(s2)), str(nl(s3)), str(nl(s4))

N_AGN_WIDE = compute_stat( t_AGN_WIDE )
N_AGN_DEEP = compute_stat( t_AGN_DEEP )
N_AGN_IR   = compute_stat( t_AGN_IR   )
N_QSO      = compute_stat( t_QSO      )
N_LyA      = compute_stat( t_LyA      )

print('any area')

print('N_AGN_WIDE', " & ", " & ".join(N_AGN_WIDE ), " \\\\")
print('N_AGN_DEEP', " & ", " & ".join(N_AGN_DEEP ), " \\\\")
print('N_AGN_IR  ', " & ", " & ".join(N_AGN_IR   ), " \\\\")
print('N_QSO     ', " & ", " & ".join(N_QSO      ), " \\\\")
print('N_LyA     ', " & ", " & ".join(N_LyA      ), " \\\\")

print('equatorial, 600 deg2')

def compute_stat_equatorial(t_survey, area=600. ):
	s0 = (t_survey['target_bit'] & 2**bitlist['AGN_WIDE'] != 0 ) & (abs(t_survey['DEC'])<5 ) & (t_survey['RA']<200 ) & (t_survey['RA']>140 )
	s1 = (t_survey['target_bit'] & 2**bitlist['AGN_DEEP'] != 0 ) & (abs(t_survey['DEC'])<5 ) & (t_survey['RA']<200 ) & (t_survey['RA']>140 )
	s2 = (t_survey['target_bit'] & 2**bitlist['AGN_IR']   != 0 ) & (abs(t_survey['DEC'])<5 ) & (t_survey['RA']<200 ) & (t_survey['RA']>140 )
	s3 = (t_survey['target_bit'] & 2**bitlist['QSO']      != 0 ) & (abs(t_survey['DEC'])<5 ) & (t_survey['RA']<200 ) & (t_survey['RA']>140 )
	s4 = (t_survey['target_bit'] & 2**bitlist['LyA']      != 0 ) & (abs(t_survey['DEC'])<5 ) & (t_survey['RA']<200 ) & (t_survey['RA']>140 )
	#print(nl(s0), nl(s1), nl(s2), nl(s3), nl(s4))
	#return str(nl(s0)), str(nl(s1)), str(nl(s2)), str(nl(s3)), str(nl(s4))
	return str(n.round(nl(s0)/area,1)), str(n.round(nl(s1)/area,1)), str(n.round(nl(s2)/area,1)), str(n.round(nl(s3)/area,1)), str(n.round(nl(s4)/area,1))

eq_N_AGN_WIDE = compute_stat_equatorial( t_AGN_WIDE )
eq_N_AGN_DEEP = compute_stat_equatorial( t_AGN_DEEP )
eq_N_AGN_IR   = compute_stat_equatorial( t_AGN_IR   )
eq_N_QSO      = compute_stat_equatorial( t_QSO      )
eq_N_LyA      = compute_stat_equatorial( t_LyA      )

print('N_AGN_WIDE', " & ", " & ".join(eq_N_AGN_WIDE ), " \\\\")
print('N_AGN_DEEP', " & ", " & ".join(eq_N_AGN_DEEP ), " \\\\")
print('N_AGN_IR  ', " & ", " & ".join(eq_N_AGN_IR   ), " \\\\")
print('N_QSO     ', " & ", " & ".join(eq_N_QSO      ), " \\\\")
print('N_LyA     ', " & ", " & ".join(eq_N_LyA      ), " \\\\")

sys.exit()

survbit_all[AGN_WIDE_all] += 2**bitlist['AGN_WIDE']
survbit_sat[AGN_WIDE_sat] += 2**bitlist['AGN_WIDE']
#
survbit_all[AGN_DEEP_all] += 2**bitlist['AGN_DEEP']
survbit_sat[AGN_DEEP_sat] += 2**bitlist['AGN_DEEP']
#
survbit_all[AGN_IR_all] += 2**bitlist['AGN_IR']
survbit_sat[AGN_IR_sat] += 2**bitlist['AGN_IR']
#
survbit_all[QSO_all] += 2**bitlist['QSO']
survbit_sat[QSO_sat] += 2**bitlist['QSO']
#
survbit_all[LyA_all] += 2**bitlist['LyA']
survbit_sat[LyA_sat] += 2**bitlist['LyA']


if os.path.isdir(dir_2_4MOST) == False:
    os.system('mkdir -p ' + dir_2_4MOST)

N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels)[::-1]:
	#HEALPIX_id = 350
	print(HEALPIX_id)
	path_2_eRO_all_catalog = os.path.join(
		dir_2_eRO_all, str(HEALPIX_id).zfill(6) + '.fit')
	path_2_eRO_sat_catalog = os.path.join(
		dir_2_eRO_sat, str(HEALPIX_id).zfill(6) + '.fit')
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
	hd_sat = fits.open(path_2_eRO_sat_catalog)
	N_agn_all = len(hd_all[1].data['ra'])
	N_agn_sat = len(hd_sat[1].data['ra'])
	#N_sat = int(N_agn_all * f_sat)
	#N_all = int(N_agn_all * (1 - f_sat))

	#rd_all = n.random.rand(N_agn_all)
	#rd_sat = n.random.rand(N_agn_sat)

	# =============================
	# =============================
	# =============================
	# eROSITA SAMPLE
	# apply flux limit erosita
	# =============================
	# =============================

	log10_fluxes = n.arange(-18., -11., 0.05)
	def detection_efficiency(log10_flux, log10_texp): return 0.5 + 0.5 * erf(
		((log10_flux + (14.0 + (log10_texp - 3.1828) * 0.5936))) / (0.3204))

	NSIDE = 512
	# healpy.nside2npix(512)
	FX_LIM = fits.open(
		os.path.join(
			os.environ['GIT_AGN_MOCK'],
			'data/erosita/flux_limits.fits'))[1].data['flux_limit_SNR3']

	# flux limit value
	hp512_all = healpy.ang2pix(
		NSIDE,
		hd_all[1].data['ra'],
		hd_all[1].data['dec'],
		nest=True,
		lonlat=True)
	FX_LIM_value_all = 10**FX_LIM[hp512_all]*1./1.

	hp512_sat = healpy.ang2pix(
		NSIDE,
		hd_sat[1].data['ra'],
		hd_sat[1].data['dec'],
		nest=True,
		lonlat=True)
	FX_LIM_value_sat = 10**FX_LIM[hp512_sat]*1./1.

	# X-ray
	detected_all = (hd_all[1].data['AGN_FX_soft'] > FX_LIM_value_all)
	detected_sat = (hd_sat[1].data['AGN_FX_soft'] > FX_LIM_value_sat)
	# optical cut
	opt_all = ( hd_all[1].data['AGN_SDSS_r_magnitude'] < 23.5 )
	opt_sat = ( hd_sat[1].data['AGN_SDSS_r_magnitude'] < 23.5 )
	# overall erosita selections
	ero_all = ( detected_all ) & ( opt_all )
	ero_sat = ( detected_sat ) & ( opt_sat )
	# optical t1
	opt_t11_all = (opt_all) & (hd_all[1].data['AGN_type']==11)
	opt_t11_sat = (opt_sat) & (hd_sat[1].data['AGN_type']==11)
	# type 2 AGN
	T2_all = (hd_all[1].data['AGN_type'] == 22) | (hd_all[1].data['AGN_type'] == 21)
	T2_sat = ( hd_sat[1].data['AGN_type'] == 22) | ( hd_sat[1].data['AGN_type'] == 21)
	# wise selection
	wise_def_all = (hd_all[1].data['WISE-W1'] < 26) & (hd_all[1].data['WISE-W1_err'] < 1) & (hd_all[1].data['WISE-W1'] < 18.85)
	wise_def_sat = (hd_sat[1].data['WISE-W1'] < 26) & (hd_sat[1].data['WISE-W1_err'] < 1) & (hd_sat[1].data['WISE-W1'] < 18.85)
	# IR sample, not detected in X-ray i.e. additional targets
	agn_ir_all =  ( detected_all == False) & (T2_all) & (opt_all) & (wise_def_all)
	agn_ir_sat =  ( detected_sat == False) & (T2_sat) & (opt_sat) & (wise_def_sat)
	# area selection
	area_erosita_all = (abs(hd_all[1].data['g_lat']) > 10) & ( hd_all[1].data['g_lon'] > 180) #& (hd_all[1].data['dec'] < 10)
	area_erosita_sat = (abs(hd_sat[1].data['g_lat']) > 10) & ( hd_sat[1].data['g_lon'] > 180) #& (hd_sat[1].data['dec'] < 10)
	# area opt-IR
	area_optIR_all = (abs(hd_all[1].data['g_lat']) > 10) #& (hd_all[1].data['dec'] < 10)
	area_optIR_sat = (abs(hd_sat[1].data['g_lat']) > 10) #& (hd_sat[1].data['dec'] < 10)
	# selection booleans
	is_erosita_all = ( (ero_all) & (area_erosita_all) )
	is_erosita_sat = ( (ero_sat) & (area_erosita_sat) )
	is_optIR_all   = ( (agn_ir_all) & (area_optIR_all) ) 
	is_optIR_sat   = ( (agn_ir_sat) & (area_optIR_sat) ) 
	is_optT11_all  = ( (opt_t11_all) & (area_optIR_all) )
	is_optT11_sat  = ( (opt_t11_sat) & (area_optIR_sat) )
	# complete selection for all 4most sub surveys
	sel_all = ( is_erosita_all ) | ( is_optIR_all ) | ( is_optT11_all )
	sel_sat = ( is_erosita_sat ) | ( is_optIR_sat ) | ( is_optT11_sat )

	# number per square degrees
	area_per_cat = healpy.nside2pixarea(8, degrees=True)
	#print(n.round(len(hd_all[1].data['redshift_R'][detected_all]) / area_per_cat,1), 'cen per deg2, Xray eROSITA detection')
	#print(n.round(len(hd_sat[1].data['redshift_R'][detected_sat]) / area_per_cat,1), 'sat per deg2, Xray eROSITA detection')
	#print('------------------------------------------------')
	##     
	#print(n.round(len(hd_all[1].data['redshift_R'][opt_all]) / area_per_cat,1), 'cen per deg2, opt r<23.5')
	#print(n.round(len(hd_sat[1].data['redshift_R'][opt_sat]) / area_per_cat,1), 'sat per deg2, opt r<23.5')
	#print('------------------------------------------------')
	##     n.round(
	#print(n.round(len(hd_all[1].data['redshift_R'][ero_all]) / area_per_cat,1), 'cen per deg2, opt + Xray')
	#print(n.round(len(hd_sat[1].data['redshift_R'][ero_sat]) / area_per_cat,1), 'sat per deg2, opt + Xray')
	#print('------------------------------------------------')
	##     n.round(
	#print(n.round(len(hd_all[1].data['redshift_R'][opt_t11_all]) / area_per_cat,1), 'cen per deg2, opt type 11')
	#print(n.round(len(hd_sat[1].data['redshift_R'][opt_t11_sat]) / area_per_cat,1), 'sat per deg2, opt type 11')
	#print('------------------------------------------------')
	##     n.round(
	#print(n.round(len(hd_all[1].data['AGN_SDSS_r_magnitude'][agn_ir_all]) / area_per_cat,1), 'cen per deg2, IR, not (opt+Xray)')
	#print(n.round(len(hd_sat[1].data['AGN_SDSS_r_magnitude'][agn_ir_sat]) / area_per_cat,1), 'sat per deg2, IR, not (opt+Xray)')
	#print('------------------------------------------------')
	##     n.round(
	#print(n.round(len(hd_all[1].data['AGN_SDSS_r_magnitude'][sel_all]) / area_per_cat,1), 'cen per deg2, all')
	#print(n.round(len(hd_sat[1].data['AGN_SDSS_r_magnitude'][sel_sat]) / area_per_cat,1), 'sat per deg2, all')
	#print('------------------------------------------------')
	# =============================
	# Define sub surveys
	# =============================
	# initialize target bits

	survbit_all = n.zeros(len(hd_all[1].data['dec']),dtype='int')
	survbit_sat = n.zeros(len(hd_sat[1].data['dec']),dtype='int')
	# subsurvey selections 
	# S6, AGN_WIDE
	AGN_WIDE_all = ( is_erosita_all ) & ( hd_all[1].data['ecl_lat'] > -80 )
	AGN_WIDE_sat = ( is_erosita_sat ) & ( hd_sat[1].data['ecl_lat'] > -80 )
	# S6, AGN_DEEP
	AGN_DEEP_all = ( is_erosita_all ) & ( hd_all[1].data['ecl_lat'] < -80 )
	AGN_DEEP_sat = ( is_erosita_sat ) & ( hd_sat[1].data['ecl_lat'] < -80 )
	# S6, AGN_IR 
	AGN_IR_all = ( is_optIR_all ) 
	AGN_IR_sat = ( is_optIR_sat ) 
	# S8, QSO
	QSO_all = ( is_optT11_all ) & ( hd_all[1].data['redshift_R'] < 2.2 )
	QSO_sat = ( is_optT11_sat ) & ( hd_sat[1].data['redshift_R'] < 2.2 )
	# S8, Lya
	LyA_all = ( is_optT11_all ) & ( hd_all[1].data['redshift_R'] > 2.2 ) 
	LyA_sat = ( is_optT11_sat ) & ( hd_sat[1].data['redshift_R'] > 2.2 ) 
	# assign target bits
	survbit_all[AGN_WIDE_all] += 2**bitlist['AGN_WIDE']
	survbit_sat[AGN_WIDE_sat] += 2**bitlist['AGN_WIDE']
	#
	survbit_all[AGN_DEEP_all] += 2**bitlist['AGN_DEEP']
	survbit_sat[AGN_DEEP_sat] += 2**bitlist['AGN_DEEP']
	#
	survbit_all[AGN_IR_all] += 2**bitlist['AGN_IR']
	survbit_sat[AGN_IR_sat] += 2**bitlist['AGN_IR']
	#
	survbit_all[QSO_all] += 2**bitlist['QSO']
	survbit_sat[QSO_sat] += 2**bitlist['QSO']
	#
	survbit_all[LyA_all] += 2**bitlist['LyA']
	survbit_sat[LyA_sat] += 2**bitlist['LyA']
	#

	target_bit = n.hstack((survbit_all[sel_all], survbit_sat[sel_sat]))


	# Now create a 4MOST file per sub survey
	# merge cen + sat
	def get_table(subsurvey):
		t_survey = Table()
		bit_value = bitlist[subsurvey]
		sel = ( target_bit & 2**bit_value != 0)

		def create_column(col_name):
			return n.hstack((hd_all[1].data[col_name][sel_all], hd_sat[1].data[col_name][sel_sat]))[sel]

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


sys.exit()

_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')
_array = create_column('')




name_array = n.array(
	[(HEALPIX_id * 1e8 + ii).astype('int').astype('str').zfill(12) for ii in n.arange(N_obj)])

t['NAME'] = Column(name_array, dtype=str)
t['RA'] = Column(ra_array, unit='degree', dtype=n.float64)
t['DEC'] = Column(dec_array, unit='degree', dtype=n.float64)
t['PMRA'] = Column(n.zeros(N_obj), unit='mas/yr', dtype=n.float32)
t['PMDEC'] = Column(n.zeros(N_obj), unit='mas/yr', dtype=n.float32)
t['EPOCH'] = Column(2000 * n.ones(N_obj), unit='yr', dtype=n.float32)
t['RESOLUTION'] = Column(n.ones(N_obj), unit='', dtype=n.int16)

s1_name = 'AGN_WIDE'
s2_name = 'AGN_DEEP'
s3_name = 'AGN_IR'


subsurvey_name = n.zeros_like(name_array).astype('str')
subsurvey_name[:] = s1_name
subsurvey_name[(ecl_lat_array < -80)] = s2_name
subsurvey_name[target_bit & 1] = s3_name

t['SUBSURVEY'] = Column(subsurvey_name, unit='', dtype=str)

priority = 100 * n.ones(N_obj)
priority[target_bit == 2] = 90
t['PRIORITY'] = Column(priority, unit='', dtype=n.int16)

template_names = n.zeros_like(z_array).astype('U100')

z_all = n.hstack((0., n.arange(0.3, 3., 0.2), 3.5, 4.5, 6.))
#z_all = n.hstack(( 0.0, n.arange(0.3, 3., 0.2), 6. ))
zmins = z_all[:-1]
zmaxs = z_all[1:]
print('ebv array', galactic_ebv_array)
ebv_1000 = (galactic_ebv_array * 1000).astype('int')
#print('ebv_galaxy', n.min(ebv_1000), n.max(ebv_1000))
ebv_1_0 = (ebv_1000 > 1000)
ebv_0_5 = (ebv_1000 > 500) & (ebv_1000 <= 1000)
ebv_0_4 = (ebv_1000 > 400) & (ebv_1000 <= 500)
ebv_0_3 = (ebv_1000 > 300) & (ebv_1000 <= 400)
ebv_0_2 = (ebv_1000 > 200) & (ebv_1000 <= 300)
ebv_0_1 = (ebv_1000 > 100) & (ebv_1000 <= 200)
ebv_0_0 = (ebv_1000 <= 100)

def z_name(z0, z1): return "_zmin_" + str(int(10 * z0)
											).zfill(2) + "_zmax_" + str(int(10 * z1)).zfill(2)

QSO = (AGN_type_array == 11) | (AGN_type_array == 12)
T2 = (AGN_type_array == 22) | (AGN_type_array == 21)
ELL = (T2) & (AGN_random_number_array < 0.2)

for z0, z1 in zip(zmins, zmaxs):
	zsel = (z_array >= z0) & (z_array < z1)
	template_names[(zsel)] = "4most_" + 'qso_BL' + \
		z_name(z0, z1) + '_EBV_0_01.fits'

	template_names[(QSO) & (zsel) & (ebv_0_0)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_01.fits'
	template_names[(QSO) & (zsel) & (ebv_0_1)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_1.fits'
	template_names[(QSO) & (zsel) & (ebv_0_2)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_2.fits'
	template_names[(QSO) & (zsel) & (ebv_0_3)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_3.fits'
	template_names[(QSO) & (zsel) & (ebv_0_4)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_4.fits'
	template_names[(QSO) & (zsel) & (ebv_0_5)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_5.fits'
	template_names[(QSO) & (zsel) & (ebv_1_0)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_1_0.fits'
	if z1 < 2.2:
		template_names[(T2) & (zsel) & (ebv_0_0)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
		template_names[(T2) & (zsel) & (ebv_0_1)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
		template_names[(T2) & (zsel) & (ebv_0_2)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
		template_names[(T2) & (zsel) & (ebv_0_3)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
		template_names[(T2) & (zsel) & (ebv_0_4)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
		template_names[(T2) & (zsel) & (ebv_0_5)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
		template_names[(T2) & (zsel) & (ebv_1_0)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'

		template_names[(ELL) & (zsel) & (ebv_0_0)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_0_01.fits'
		template_names[(ELL) & (zsel) & (ebv_0_1)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_0_1.fits'
		template_names[(ELL) & (zsel) & (ebv_0_2)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_0_2.fits'
		template_names[(ELL) & (zsel) & (ebv_0_3)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_0_3.fits'
		template_names[(ELL) & (zsel) & (ebv_0_4)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_0_4.fits'
		template_names[(ELL) & (zsel) & (ebv_0_5)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_0_5.fits'
		template_names[(ELL) & (zsel) & (ebv_1_0)] = "4most_" + \
			'LRG' + z_name(z0, z1) + '_EBV_1_0.fits'
	if z1 >= 2.2 and z1 < 6.:
		template_names[(T2) & (zsel) & (ebv_0_0)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
		template_names[(T2) & (zsel) & (ebv_0_1)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
		template_names[(T2) & (zsel) & (ebv_0_2)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
		template_names[(T2) & (zsel) & (ebv_0_3)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
		template_names[(T2) & (zsel) & (ebv_0_4)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
		template_names[(T2) & (zsel) & (ebv_0_5)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
		template_names[(T2) & (zsel) & (ebv_1_0)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'

		template_names[(ELL) & (zsel) & (ebv_0_0)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
		template_names[(ELL) & (zsel) & (ebv_0_1)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
		template_names[(ELL) & (zsel) & (ebv_0_2)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
		template_names[(ELL) & (zsel) & (ebv_0_3)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
		template_names[(ELL) & (zsel) & (ebv_0_4)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
		template_names[(ELL) & (zsel) & (ebv_0_5)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
		template_names[(ELL) & (zsel) & (ebv_1_0)] = "4most_" + \
			'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'

zsel = (z_array > 6.)
if len(zsel.nonzero()[0]) > 0:
	#print(z0, z1)
	template_names[(QSO) & (zsel) & (ebv_0_0)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_01.fits'
	template_names[(QSO) & (zsel) & (ebv_0_1)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_1.fits'
	template_names[(QSO) & (zsel) & (ebv_0_2)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_2.fits'
	template_names[(QSO) & (zsel) & (ebv_0_3)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_3.fits'
	template_names[(QSO) & (zsel) & (ebv_0_4)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_4.fits'
	template_names[(QSO) & (zsel) & (ebv_0_5)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_0_5.fits'
	template_names[(QSO) & (zsel) & (ebv_1_0)] = "4most_" + \
		'qso_BL' + z_name(z0, z1) + '_EBV_1_0.fits'
	template_names[(T2) & (zsel) & (ebv_0_0)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
	template_names[(T2) & (zsel) & (ebv_0_1)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
	template_names[(T2) & (zsel) & (ebv_0_2)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
	template_names[(T2) & (zsel) & (ebv_0_3)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
	template_names[(T2) & (zsel) & (ebv_0_4)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
	template_names[(T2) & (zsel) & (ebv_0_5)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
	template_names[(T2) & (zsel) & (ebv_1_0)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'
	template_names[(ELL) & (zsel) & (ebv_0_0)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_01.fits'
	template_names[(ELL) & (zsel) & (ebv_0_1)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_1.fits'
	template_names[(ELL) & (zsel) & (ebv_0_2)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_2.fits'
	template_names[(ELL) & (zsel) & (ebv_0_3)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_3.fits'
	template_names[(ELL) & (zsel) & (ebv_0_4)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_4.fits'
	template_names[(ELL) & (zsel) & (ebv_0_5)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_0_5.fits'
	template_names[(ELL) & (zsel) & (ebv_1_0)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'
	template_names[(ELL) & (zsel) & (ebv_1_0)] = "4most_" + \
		'AGN_type2' + z_name(z0, z1) + '_EBV_1_0.fits'

# sanity checks

tpls = sorted(n.array(list(set(template_names))))
##print(tpls)
print('N templates used', len(tpls))
bad = (template_names == '0.0')
print(len(bad.nonzero()[0]))

N_all = len(bad)


t['TEMPLATE'] = Column(template_names, unit='', dtype=str)

ruleset_array = n.zeros_like(z_array).astype('U20')
ruleset_array[ruleset_array == "0.0"] = "AGN_ALL_3PC"
#print(template_names, ruleset_array)
#print(template_names.shape, ruleset_array.shape)
#print(template_names.dtype, ruleset_array.dtype)
t['RULESET'] = Column(ruleset_array, unit='', dtype=str)
t['REDSHIFT_ESTIMATE'] = Column(z_array, unit='', dtype=n.float32)
t['REDSHIFT_ERROR'] = Column(
	n.ones(N_obj) *
	0.00001,
	unit='',
	dtype=n.float32)

t['EXTENT_PARAMETER'] = Column(n.ones(N_obj), unit='', dtype=n.float32)
t['REDDENING'] = Column(galactic_ebv_array, unit='mag', dtype=n.float32)
t['DATE_EARLIEST'] = Column(
	n.ones(N_obj) * 59215.,
	unit='',
	dtype=n.float64)
t['DATE_LATEST'] = Column(
	n.ones(N_obj) * 66520,
	unit='',
	dtype=n.float64)

magnitude_4fs = n.hstack(
	(hd_all[1].data['AGN_SDSS_r_magnitude'][sel_all],
		hd_sat[1].data['AGN_SDSS_r_magnitude'][sel_sat]))
magnitude_4fs_err = n.ones(N_obj) * 0.1
# implement the Palanque Delabrouille fraction to mimic fibermag !
zi, zf, dm = n.loadtxt(
	os.path.join(
		os.environ['GIT_VS'], 'data/m_psf-m_model_redshift_PD16.txt'), unpack=True)
for z0i, z1i, dmi in zip(zi, zf, dm):
	#print(z0i, z1i)
	ssel = (z_array >= z0i) & (z_array < z1i)
	magnitude_4fs[ssel] = magnitude_4fs[ssel] + \
		dmi * n.ones_like(magnitude_4fs[ssel])

mag_type = n.zeros(N_obj).astype('U10')
mag_type[mag_type == "0.0"] = "SDSS_r_AB"

t['MAG'] = Column(magnitude_4fs, unit='mag', dtype=n.float32)
t['MAG_ERR'] = Column(magnitude_4fs_err, unit='', dtype=n.float32)
t['MAG_TYPE'] = Column(mag_type, unit='', dtype=n.str)

print(path_2_4MOST_catalog)
if os.path.isfile(path_2_4MOST_catalog):
	os.system("rm " + path_2_4MOST_catalog)

#t.write(path_2_4MOST_catalog, format='fits')
