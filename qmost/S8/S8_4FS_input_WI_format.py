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

survey = 'S8'
working_dir = os.path.join(os.environ['HOME'], 'data', '4most', survey)
# name of the sub survey for which the LSM is created
sub_survey_names = np.array([ 'BG', 'LRG', 'ELG', 'QSO', 'LyA'])
# loads the catalog
catalog_input = os.path.join(working_dir, 'S8_Cosmology_latest_input_cat.fits.gz')
data = fits.open(catalog_input)[1].data

path_2_out = catalog_input[:-8]+'_4FS_new_format.fits'


#======
# vertical stripe in the ngc between kid-s and atlas; to be removed for all targets
ra230stripe = (data['ra']>229) & (data['ra']<233) & (data['dec']>-10) & (data['dec']<-5)
# bg
s1 = (data['SUB_CAT_ID']==1) & (ra230stripe==False)
# lrg
ind = np.where((data['SUB_CAT_ID']==2) & (ra230stripe==False))[0]
ind = np.random.choice(ind,size=int(0.999*len(ind)),replace=False) # subsampling by 0.571, to adapt density from 700/deg2 to 400/deg2
s2 = np.zeros(len(data),dtype=bool)
s2[ind] = True
# elg [footprint = waves-wide: -30<ra<50 , -36<dec<-26 ; and ~300 deg2 ~inside SPT: 0<ra<50,-45<dec<-46]
elgfoot = ( (((data['ra']>330) | (data['ra']<50)) & (data['dec']>-36) & (data['dec']<-26)) | ((data['ra']>0) & (data['ra']<50) & (data['dec']>-45) & (data['dec']<-36)) )
ind = np.where((data['SUB_CAT_ID']==3) & (elgfoot) & (ra230stripe==False))[0]
ind = np.random.choice(ind,size=int(0.999*len(ind)),replace=False) # subsampling by 0.923, to adapt density from 1300/deg2 to 1200/deg2
s3 = np.zeros(len(data),dtype=bool)
s3[ind] = True
# qso
s4 = (data['SUB_CAT_ID']==4) & (ra230stripe==False)
# lya
s5 = (data['SUB_CAT_ID']==5) & (ra230stripe==False)
#======

N_tot = len(data['IDNUM']  )

sub_cats = np.zeros_like(data['REDSHIFT']).astype('U10')
sub_cats[s1] = sub_survey_names[0]
sub_cats[s2] = sub_survey_names[1]
sub_cats[s3] = sub_survey_names[2]
sub_cats[s4] = sub_survey_names[3]
sub_cats[s5] = sub_survey_names[4]

mag_type = np.zeros_like(data['REDSHIFT']).astype('U9')
mag_type[mag_type=="0.0"] = "SDSS_r_AB"
mag_type2 = np.zeros((1,N_tot)).astype('U9')
mag_type2[0]=mag_type

ruleset = data['RULESET']
ruleset[s5] = "COSMO_LYA"

#### NEED TO REDEFINE AREAS PROPERLY.

ok =  ((s1) | (s2) | (s3) | (s4) | (s5) ) & (data['R_MAG']>15)

from astropy.table import Table, Column
t = Table()

t['NAME'              ] = Column( np.arange(N_tot).astype('U7') [ok] )
t['RA'                ] = Column( data['RA']                     [ok] )   
t['DEC'               ] = Column( data['DEC']                    [ok] )   
t['PMRA'              ] = Column( np.zeros(N_tot)                [ok] )   
t['PMDEC'             ] = Column( np.zeros(N_tot)                [ok] )   
t['EPOCH'             ] = Column( np.ones(N_tot)          [ok] *2000.0)   
t['RESOLUTION'        ] = Column( np.ones((1,N_tot)).T           [ok] ) 
t['SUBSURVEY'         ] = Column( sub_cats                       [ok] )  
t['PRIORITY'          ] = Column( np.ones(N_tot)           [ok] *100.0)
t['TEMPLATE'          ] = Column( data['TEMPLATE']               [ok] )
t['RULESET'           ] = Column( ruleset                        [ok] ) 
t['REDSHIFT_ESTIMATE' ] = Column( data['REDSHIFT']               [ok] )
t['REDSHIFT_ERROR'    ] = Column( np.ones(N_tot)         [ok] *0.00001)
t['EXTENT_PARAMETER'  ] = Column( np.ones(N_tot)                 [ok] )   
t['MAG'               ] = Column( data['R_MAG']                  [ok] )   
t['MAG_ERR'           ] = Column( np.ones(N_tot)            [ok] *0.01)   
t['MAG_TYPE'          ] = Column( mag_type2.T                    [ok] )   
t['REDDENING'         ] = Column( data['REDSHIFT']            [ok]*0. )
t['DATE_EARLIEST'     ] = Column( np.ones(N_tot)          [ok]*59215. )  
t['DATE_LATEST'       ] = Column( np.ones(N_tot)          [ok]*66520  )   

print( path_2_out )
if os.path.isfile(path_2_out):
	os.system("rm "+path_2_out)

t.write(path_2_out, format='fits')

sys.exit()

hdu_cols = fits.ColDefs([
fits.Column(name='NAME'              , format='20A',  array=np.arange(N_tot).astype('str') ),
fits.Column(name='RA'                , format='1D' ,  array = ra [eRO_DE] ),   
fits.Column(name='DEC'               , format='1D' ,  array = dec[eRO_DE] ),   
fits.Column(name='PMRA'              , format='1E' ,  array = np.zeros(N_tot) ),   
fits.Column(name='PMDEC'             , format='1E' ,  array = np.zeros(N_tot) ),   
fits.Column(name='EPOCH'             , format='1E' ,  array = np.ones(N_tot)*2000.0 ),   
fits.Column(name='RESOLUTION'        , format='1I' ,  array = np.ones(N_tot) ), 
fits.Column(name='SUBSURVEY'         , format='10A' , array = sub_cats[eRO_DE] ),  
fits.Column(name='PRIORITY'          , format='1I' ,  array = np.ones(N_tot)*100.0 ),
fits.Column(name='TEMPLATE'          , format='100A', array=template_names[eRO_DE] ),
fits.Column(name='RULESET'           , format='20A' , array=ruleset_array[eRO_DE] ), 
fits.Column(name='REDSHIFT_ESTIMATE' , format='1E' ,  array=zs[eRO_DE] ),
fits.Column(name='REDSHIFT_ERROR'    , format='1E' ,  array = np.ones(N_tot)*0.00001 ),
fits.Column(name='EXTENT_PARAMETER'  , format='1E' ,  array = np.ones(N_tot) ),   
fits.Column(name='MAG'               , format='1E' ,  array = magnitude_4fs[eRO_DE] ),   
fits.Column(name='MAG_ERR'           , format='1E' ,  array = magnitude_4fs_err[eRO_DE] ),   
fits.Column(name='MAG_TYPE'          , format='55A',  array = mag_type2 ),   
fits.Column(name='REDDENING'         , format='1E' ,  array=data[1].data['ebv_galaxy'][eRO_DE]),
fits.Column(name='DATE_EARLIEST'     , format='1D' ,  array=np.ones(N_tot)*59215. ),   
fits.Column(name='DATE_LATEST'       , format='1D' ,  array=np.ones(N_tot)*66520 )   
])

tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
#define the header
prihdr = fits.Header()
prihdr['author'] = 'COMPARAT'
prihdu = fits.PrimaryHDU(header=prihdr)
#writes the file
thdulist = fits.HDUList([prihdu, tb_hdu])

print( path_2_out )
if os.path.isfile(path_2_out):
	os.system("rm "+path_2_out)
thdulist.writeto(path_2_out)

os.system('gzip '+path_2_out)