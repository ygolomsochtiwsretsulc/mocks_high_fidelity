"""
4MOST 4FS_WI SSM script

python script

Generate a fits file containing an example Small Scale Merit SSM

python dependencies:
 - os, sys, numpy, scipy
 - astropy

> pip install ...

It creates a SSM file that interpolates the SSM function required for a sub survey.

"""
import os, sys
import numpy as np

import astropy.io.fits as fits

survey = 'S8'
working_dir = os.path.join(os.environ['HOME'], 'data', '4most', survey)
# name of the sub survey for which the LSM is created
sub_survey_names = np.array([ 'BG', 'LRG', 'ELG', 'QSO', 'LyA'])
# loads the catalog
t_exp_max = np.array([ 40. ,
   40. ,
   30. ,
  60. ,
  60. 
 ])
goal_area =  np.array([  7500., #'COSMO_BG',
   7500., #'COSMO_LRG',
   1000., #'COSMO_ELG',
   7500., #'COSMO_QSO',
   7500. #'COSMO_LyAlpha',
  ])


#sub_survey_name = sub_survey_names[0]
#tmax_value = t_exp_max[0]
#area_required = goal_area[0]

################################################
################################################
################################################
# AREQ, TMAX
################################################
################################################
################################################

def create_files(sub_survey_name, tmax_value, area_required):
	# name of the output file
	out_file_TMAX = os.path.join(working_dir, survey +'_TMAX_'+str(sub_survey_name)+'.fits')
	out_file_AREQ = os.path.join(working_dir, survey +'_AREQ_'+str(sub_survey_name)+'.fits')
	#plot_file = os.path.join(working_dir,survey +'_LSM_'+str(sub_survey_name)+'.png')

	cols = fits.ColDefs([
		fits.Column( "SUBSURVEY"	 ,unit='' ,format="20A", array=np.array([sub_survey_name])), 
		fits.Column( "TMAX"            ,unit='' ,format="D", array=np.array([tmax_value]) ) 
		])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.header['author'] = 'JC'
	tbhdu.header['HIERARCH SUBSURVEY'] = sub_survey_name
	print( out_file_TMAX )
	if os.path.isfile(out_file_TMAX):
		os.system("rm "+out_file_TMAX)
	tbhdu.writeto(out_file_TMAX)

	cols = fits.ColDefs([
		fits.Column( "SUBSURVEY"	 ,unit='' ,format="20A", array=np.array([sub_survey_name])), 
		fits.Column( "AREA_REQ"            ,unit='' ,format="D", array=np.array([area_required]) ) 
		])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.header['author'] = 'JC'
	tbhdu.header['HIERARCH SUBSURVEY'] = sub_survey_name
	print( out_file_AREQ )
	if os.path.isfile(out_file_AREQ):
		os.system("rm "+out_file_AREQ)
	tbhdu.writeto(out_file_AREQ)

for sub_survey_name, tmax_value, area_required in zip(sub_survey_names,t_exp_max, goal_area ):
	create_files(sub_survey_name, tmax_value, area_required)

lsm_files  = np.array([ survey +'_LSM_'+str(sub_survey_name)+'.fits' for sub_survey_name in sub_survey_names ])
ssm_files  = np.array([ survey +'_SSM_'+str(sub_survey_name)+'.fits' for sub_survey_name in sub_survey_names ])
tmax_files = np.array([ survey +'_TMAX_'+str(sub_survey_name)+'.fits' for sub_survey_name in sub_survey_names ])
areq_files = np.array([ survey +'_AREQ_'+str(sub_survey_name)+'.fits' for sub_survey_name in sub_survey_names ])

out_file_PARA = os.path.join(working_dir, survey +'_SUBSURVEY_PARAMS.fits')

cols = fits.ColDefs([
	fits.Column( "SUBSURVEY"	 ,unit='' ,format="20A", array = sub_survey_names), 
	fits.Column( "LSM_FILENAME"  ,unit='' ,format="200A", array  = lsm_files  ), 
	fits.Column( "SSM_FILENAME"  ,unit='' ,format="200A", array  = ssm_files  ), 
	fits.Column( "TMAX_FILENAME"  ,unit='' ,format="200A", array = tmax_files ), 
	fits.Column( "AREQ_FILENAME"  ,unit='' ,format="200A", array = areq_files ), 
	])
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.header['author'] = 'JC'
print( out_file_PARA )
if os.path.isfile(out_file_PARA):
	os.system("rm "+out_file_PARA)
tbhdu.writeto(out_file_PARA)
 

