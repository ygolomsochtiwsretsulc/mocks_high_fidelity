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

env = sys.argv[1] # 'MD10'
root_dir = os.path.join(os.environ[env])

nl = lambda selection : len(selection.nonzero()[0])

dir_2_AGN_all = os.path.join(root_dir, "AGN_all.fits")
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
