"""
S5 SSM

It creates a SSM file that interpolates the SSM function required for a sub survey.

"""
import os, sys
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
import astropy.io.fits as fits

survey = 'S8'
working_dir = os.path.join(os.environ['HOME'], 'data', '4most', survey)
# name of the sub survey for which the LSM is created
sub_survey_names = np.array([ 'BG', 'LRG', 'ELG', 'QSO', 'LyA'])

def create_SSM_file(sub_survey_name):
	# name of the output file
	out_file = os.path.join(working_dir, survey +'_SSM_'+str(sub_survey_name)+'.fits')
	plot_file = os.path.join(working_dir,survey +'_SSM_'+str(sub_survey_name)+'.png')
	#########
	#########
	# 1. First create the grid of pixels to interpolate on
	#########
	#########
	x = np.arange(-0.01, 1.02, 0.01)
	ssm_values = np.zeros(len(x))
	#########
	#########
	# 2. Choose the SSM function that fits your needs
	#########
	#########
	if sub_survey_name == 'BG' : 
		ssm_fun = lambda x: 0.5+0.5*erf((x-0.8)/0.1)

	if sub_survey_name == 'LRG' : 
		ssm_fun = lambda x: 0.5+0.5*erf((x-0.8)/0.1)

	if sub_survey_name == 'ELG' : 
		ssm_fun = lambda x: 0.5+0.5*erf((x-0.8)/0.1) # x #**(1/2.) # 0.5+0.5*erf((x-0.8)/0.05)

	if sub_survey_name == 'QSO' : 
		ssm_fun = lambda x: 0.5+0.5*erf((x-0.8)/0.1) # x #**(1/2.) # x**7.

	if sub_survey_name == 'LyA' : 
		ssm_fun = lambda x: 0.5+0.5*erf((x-0.8)/0.1) # x #**(1/2.) # x**7.

	#########
	#########
	# 2. Assigns SSM values
	#########
	#########
	ssm_values = ssm_fun(x)
	print('SSM values:',ssm_fun(np.arange(0.,1.1,0.1)))
	cols = fits.ColDefs([
		fits.Column( "completeness"	 ,unit='' ,format="D", array=x), 
		fits.Column( "SSM"            ,unit='' ,format="D", array=ssm_values ) 
		])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.header['author'] = 'JC'
	tbhdu.header['HIERARCH SUBSURVEY'] = sub_survey_name
	if os.path.isfile(out_file):
		os.remove(out_file)
	tbhdu.writeto(out_file)
	# Creates the Figure
	plt.figure(0, (6,6))
	plt.plot(x, ssm_values)
	plt.grid()
	plt.ylim((-0.05,1.05))
	plt.xlim((-0.05,1.05))
	plt.xlabel('Completeness (fraction)')
	plt.ylabel('Small scale merit value')
	plt.title(sub_survey_name)
	plt.savefig(plot_file)
	plt.clf()

for sub_survey_name in sub_survey_names:
	print(sub_survey_name)
	create_SSM_file(sub_survey_name)