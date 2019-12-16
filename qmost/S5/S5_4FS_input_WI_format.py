"""
Creates catalog compliant with the web interface
Creates a fits catalog containing the 4FS input columns.

Create the ID array
remove unwanted column
check rulesets and templates names are correct 

"""
print('CREATES 4FS WI FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')

import os
import sys
import astropy.io.fits as fits
import numpy as np
from astropy.coordinates import SkyCoord

#from lib_model_agn import *

survey = 'S5'
working_dir = os.path.join(os.environ['HOME'], 'data', '4most', survey)
# name of the sub survey for which the LSM is created
sub_survey_names = np.array([ 'cluster_BCG', 'cluster_redGAL', 'filament_GAL'])
# loads the catalog
catalog_input = os.path.join(working_dir, 'S5_Clusters_latest_input_cat.fits.gz')
data = fits.open(catalog_input)[1].data

sub_survey_names = np.array([ 'cluster_BCG', 'cluster_redGAL', 'filament_GAL'])

s1 = (data['ID_SUB_SURVEY']==1)
s2 = (data['ID_SUB_SURVEY']==2)
s3 = (data['ID_SUB_SURVEY']==3)

path_2_out = catalog_input[:-8]+'_4FS_new_format.fits'

N_tot = len(data['IDNUM']  )

sub_cats = np.zeros_like(data['REDSHIFT']).astype('U20')
sub_cats[s1] = sub_survey_names[0]
sub_cats[s2] = sub_survey_names[1]
sub_cats[s3] = sub_survey_names[2]

mag_type = np.zeros_like(data['REDSHIFT']).astype('U9')
mag_type[mag_type=="0.0"] = "SDSS_r_AB"
mag_type2 = np.zeros((1,N_tot)).astype('U9')
mag_type2[0]=mag_type

ok =(data['DEC']<5)

from astropy.table import Table, Column
t = Table()

t['NAME'              ] = Column( np.arange(N_tot).astype('str')[ok] )
t['RA'                ] = Column( data['RA'][ok] )   
t['DEC'               ] = Column( data['DEC'][ok] )   
t['PMRA'              ] = Column( np.zeros(N_tot)[ok] )   
t['PMDEC'             ] = Column( np.zeros(N_tot)[ok] )   
t['EPOCH'             ] = Column( np.ones(N_tot)[ok]*2000.0 )   
t['RESOLUTION'        ] = Column( np.ones((1,N_tot)).T[ok] ) 
t['SUBSURVEY'         ] = Column( sub_cats[ok] )  
t['PRIORITY'          ] = Column( np.ones(N_tot)[ok]*100.0 )
t['TEMPLATE'          ] = Column( data['TEMPLATE'][ok] )
t['RULESET'           ] = Column( data['RULESET'][ok] ) 
t['REDSHIFT_ESTIMATE' ] = Column( data['REDSHIFT'][ok] )
t['REDSHIFT_ERROR'    ] = Column( np.ones(N_tot)[ok]*0.00001 )
t['EXTENT_PARAMETER'  ] = Column( np.ones(N_tot)[ok] )   
t['MAG'               ] = Column( data['R_MAG'][ok] )   
t['MAG_ERR'           ] = Column( np.ones(N_tot)[ok]*0.01   )   
t['MAG_TYPE'          ] = Column( mag_type2.T[ok]           )   
t['REDDENING'         ] = Column( data['EBV'][ok]           )
t['DATE_EARLIEST'     ] = Column( np.ones(N_tot)[ok]*59215. )  
t['DATE_LATEST'       ] = Column( np.ones(N_tot)[ok]*66520  )   

print( path_2_out )
if os.path.isfile(path_2_out):
	os.system("rm "+path_2_out)

t.write(path_2_out, format='fits')
