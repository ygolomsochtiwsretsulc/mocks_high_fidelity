"""
Creates a fits catalog containing the 4FS input columns.

4FS input catalog columns

     1  OBJECT_ID                  15A
     2  IDNUM                      1J
     3  RA                         1D    
     4  DEC                        1D    
     5  PRIORITY                   1I
     6  ORIG_TEXP_B                1E    
     7  ORIG_TEXP_D                1E    
     8  ORIG_TEXP_G                1E    
     9  RESOLUTION                 1B
     10 R_MAG                     1E    
     11 TEMPLATE                  30A
     12 RULESET                   9A

Create the ID array
remove unwanted column
check rulesets and templates names are correct 

"""
print('CREATES 4FS FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')

import os
import sys
import astropy.io.fits as fits
import numpy as n
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column

survey = 'S6'
working_dir = os.path.join(os.environ['HOME'], 'data', '4most', survey)
catalog_input = os.path.join(working_dir, 'eRosita_SNR3_with_photometry.fits.gz')
data = fits.open(catalog_input)

catalog_input_3 = os.path.join(working_dir, 'all_IR_AGN.fits.gz')
data_3 = fits.open(catalog_input_3)

path_2_out = 'eROSITA_component_4FS_new_format.fits'
path_2_out_3 = 'IR_component_4FS_new_format.fits'

zs = data[1].data['redshift_R']
print('min redshift ',n.min(zs))

mag = 'SDSS_r_magnitude'
magnitude_4fs=data[1].data[mag]
magnitude_4fs_err=data[1].data['HSC-r_err']

template_names = n.zeros_like(data[1].data['redshift_R']).astype('U100')
ruleset_array = n.zeros_like(data[1].data['redshift_R']).astype('U20')
ruleset_array[ruleset_array=="0.0"] = "AGN_ALL_3PC"
print(ruleset_array) 


# reassigns templates correctly
z_all = n.hstack(( 0., n.arange(0.3, 3., 0.2), 3.5, 4.5, 6. ))
#z_all = n.hstack(( 0.0, n.arange(0.3, 3., 0.2), 6. ))
zmins = z_all[:-1]
zmaxs = z_all[1:]

ebv_1000 = (data[1].data['ebv_galaxy']*1000).astype('int')
print('ebv_galaxy', n.min(ebv_1000), n.max(ebv_1000))
ebv_1_0 = ( ebv_1000 > 1000 ) 
ebv_0_5 = ( ebv_1000 > 500 ) & ( ebv_1000 <= 1000 ) 
ebv_0_4 = ( ebv_1000 > 400 ) & ( ebv_1000 <= 500 ) 
ebv_0_3 = ( ebv_1000 > 300 ) & ( ebv_1000 <= 400 ) 
ebv_0_2 = ( ebv_1000 > 200 ) & ( ebv_1000 <= 300 ) 
ebv_0_1 = ( ebv_1000 > 100 ) & ( ebv_1000 <= 200 ) 
ebv_0_0 = ( ebv_1000 <= 100 ) 

z_name = lambda z0, z1 : "_zmin_"+str(int(10*z0)).zfill(2)+"_zmax_"+str(int(10*z1)).zfill(2)

agn_type = data[1].data['AGN_type']
QSO = (agn_type==11)|(agn_type==12)
T2  = (agn_type==22)|(agn_type==21)
ELL = (T2) & (data[1].data['AGN_random_number']<0.2)
for z0,z1 in zip(zmins,zmaxs):
	zsel = (zs>=z0)&(zs<z1)
	template_names[(zsel)]                     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_01.fits'  
	
	template_names[(QSO)&(zsel)&(ebv_0_0)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_01.fits'  
	template_names[(QSO)&(zsel)&(ebv_0_1)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_1.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_2)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_2.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_3)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_3.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_4)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_4.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_5)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_5.fits'   
	template_names[(QSO)&(zsel)&(ebv_1_0)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_1_0.fits'   
	if z1<2.2:
		template_names[(T2)&(zsel)&(ebv_0_0)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_01.fits'
		template_names[(T2)&(zsel)&(ebv_0_1)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_1.fits' 
		template_names[(T2)&(zsel)&(ebv_0_2)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_2.fits' 
		template_names[(T2)&(zsel)&(ebv_0_3)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_3.fits' 
		template_names[(T2)&(zsel)&(ebv_0_4)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_4.fits' 
		template_names[(T2)&(zsel)&(ebv_0_5)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_5.fits' 
		template_names[(T2)&(zsel)&(ebv_1_0)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_1_0.fits' 
		
		template_names[(ELL)&(zsel)&(ebv_0_0)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_01.fits'  
		template_names[(ELL)&(zsel)&(ebv_0_1)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_1.fits'   
		template_names[(ELL)&(zsel)&(ebv_0_2)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_2.fits'   
		template_names[(ELL)&(zsel)&(ebv_0_3)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_3.fits'   
		template_names[(ELL)&(zsel)&(ebv_0_4)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_4.fits'   
		template_names[(ELL)&(zsel)&(ebv_0_5)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_5.fits'   
		template_names[(ELL)&(zsel)&(ebv_1_0)]    = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_1_0.fits'   
	if z1>=2.2 and z1<6.:
		template_names[(T2)&(zsel)&(ebv_0_0)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_01.fits'
		template_names[(T2)&(zsel)&(ebv_0_1)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_1.fits' 
		template_names[(T2)&(zsel)&(ebv_0_2)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_2.fits' 
		template_names[(T2)&(zsel)&(ebv_0_3)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_3.fits' 
		template_names[(T2)&(zsel)&(ebv_0_4)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_4.fits' 
		template_names[(T2)&(zsel)&(ebv_0_5)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_5.fits' 
		template_names[(T2)&(zsel)&(ebv_1_0)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_1_0.fits' 
		
		template_names[(ELL)&(zsel)&(ebv_0_0)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_01.fits'
		template_names[(ELL)&(zsel)&(ebv_0_1)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_1.fits' 
		template_names[(ELL)&(zsel)&(ebv_0_2)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_2.fits' 
		template_names[(ELL)&(zsel)&(ebv_0_3)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_3.fits' 
		template_names[(ELL)&(zsel)&(ebv_0_4)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_4.fits' 
		template_names[(ELL)&(zsel)&(ebv_0_5)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_5.fits' 
		template_names[(ELL)&(zsel)&(ebv_1_0)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_1_0.fits' 

zsel = (zs>6.)
if len(zsel.nonzero()[0])>0:
	print(z0,z1)
	template_names[(QSO)&(zsel)&(ebv_0_0)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_01.fits'  
	template_names[(QSO)&(zsel)&(ebv_0_1)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_1.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_2)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_2.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_3)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_3.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_4)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_4.fits'   
	template_names[(QSO)&(zsel)&(ebv_0_5)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_0_5.fits'   
	template_names[(QSO)&(zsel)&(ebv_1_0)]     = "4most_"+'qso_BL'+z_name( z0, z1)+'_EBV_1_0.fits'   
	template_names[(T2)&(zsel)&(ebv_0_0)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_01.fits'
	template_names[(T2)&(zsel)&(ebv_0_1)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_1.fits' 
	template_names[(T2)&(zsel)&(ebv_0_2)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_2.fits' 
	template_names[(T2)&(zsel)&(ebv_0_3)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_3.fits' 
	template_names[(T2)&(zsel)&(ebv_0_4)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_4.fits' 
	template_names[(T2)&(zsel)&(ebv_0_5)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_5.fits' 
	template_names[(T2)&(zsel)&(ebv_1_0)]     = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_1_0.fits' 
	template_names[(ELL)&(zsel)&(ebv_0_0)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_01.fits'
	template_names[(ELL)&(zsel)&(ebv_0_1)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_1.fits' 
	template_names[(ELL)&(zsel)&(ebv_0_2)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_2.fits' 
	template_names[(ELL)&(zsel)&(ebv_0_3)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_3.fits' 
	template_names[(ELL)&(zsel)&(ebv_0_4)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_4.fits' 
	template_names[(ELL)&(zsel)&(ebv_0_5)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_0_5.fits' 
	template_names[(ELL)&(zsel)&(ebv_1_0)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_1_0.fits' 
	template_names[(ELL)&(zsel)&(ebv_1_0)]    = "4most_"+'AGN_type2'+z_name( z0, z1)+'_EBV_1_0.fits' 

# sanity checks

tpls = n.array(list(set(template_names)))
tpls.sort()
print(tpls)
print('N templates used', len(tpls))
bad = (template_names=='0.0')
print(len(bad.nonzero()[0]))

N_all = len(bad)

print(template_names, ruleset_array)
print(template_names.shape, ruleset_array.shape)
print(template_names.dtype, ruleset_array.dtype)



# implement the Palanque Delabrouille fraction to mimic fibermag !
zi,zf,dm=n.loadtxt(os.path.join(os.environ['GIT_VS'],'data/m_psf-m_model_redshift_PD16.txt'), unpack=True)
for z0i, z1i, dmi in zip(zi,zf,dm):
	ssel = (zs>=z0i)&(zs<z1i)
	magnitude_4fs[ssel] = magnitude_4fs[ssel] + dmi*n.ones_like(magnitude_4fs[ssel])

ra = data[1].data['ra']
dec = data[1].data['dec']
coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
bb = coords.galactic.b.value
ll = coords.galactic.l.value
bb_ecl = coords.barycentrictrueecliptic.lat
ll_ecl = coords.barycentrictrueecliptic.lon
sel = (bb_ecl.value<-80)

area_sel = (ll>180) & (abs(bb)>20) & (dec<5) & (dec>-80) 
s1 = (area_sel) & (sel==False) & (magnitude_4fs < 23.0) # 22.8
s2 = (area_sel) & (sel) & (magnitude_4fs < 23.5) # 23.2
eRO_DE =  (s1) | (s2) 

#AGN_ALL_*PC_rXXY = ruleset to be used for all targets with r.gt.22.0, with magnitude XX.Y<r<=XX.Y+0.1
#faint_ruleset = n.array(['AGN_ALL_3PC_r'+str(int(rmag*10)) for rmag in magnitude_4fs ])
#ruleset_array[magnitude_4fs>22] = faint_ruleset[magnitude_4fs>22]

s1_name = 'AGN_WIDE'
s2_name = 'AGN_DEEP'

sub_cats = n.zeros_like(data[1].data['redshift_R']).astype('U20')
sub_cats[s1] = s1_name
sub_cats[s2] = s2_name


N_tot = len(ra[eRO_DE])

mag_type = n.zeros_like(data[1].data['redshift_R']).astype('U9')
mag_type[mag_type=="0.0"] = "SDSS_r_AB"
mag_type2 = n.zeros((1,len(data[1].data['redshift_R'][eRO_DE]))).astype('U9')
mag_type2[0]=mag_type[eRO_DE]
N_test = len(mag_type2[0])

magnitude_4fs2 = n.zeros((1,len(data[1].data['redshift_R'][eRO_DE])))
magnitude_4fs2[0] = magnitude_4fs[eRO_DE]

magnitude_4fs_err2 = n.zeros((1,len(data[1].data['redshift_R'][eRO_DE])))
magnitude_4fs_err2[0] = magnitude_4fs_err[eRO_DE]


t = Table()

t['NAME'              ] = Column( n.arange(N_tot).astype('str') )
t['RA'                ] = Column( ra [eRO_DE] )   
t['DEC'               ] = Column( dec[eRO_DE] )   
t['PMRA'              ] = Column( n.zeros(N_tot) )   
t['PMDEC'             ] = Column( n.zeros(N_tot) )   
t['EPOCH'             ] = Column( n.ones(N_tot)*2000.0 )   
t['RESOLUTION'        ] = Column( n.ones((1,N_tot)).T ) 
t['SUBSURVEY'         ] = Column( sub_cats[eRO_DE] )  
t['PRIORITY'          ] = Column( n.ones(N_tot)*100.0 )
t['TEMPLATE'          ] = Column( template_names[eRO_DE] )
t['RULESET'           ] = Column( ruleset_array[eRO_DE] ) 
t['REDSHIFT_ESTIMATE' ] = Column( zs[eRO_DE] )
t['REDSHIFT_ERROR'    ] = Column( n.ones(N_tot)*0.00001 )
t['EXTENT_PARAMETER'  ] = Column( n.ones(N_tot) )   
t['MAG'               ] = Column( magnitude_4fs2.T )   
t['MAG_ERR'           ] = Column( magnitude_4fs_err2.T )   
t['MAG_TYPE'          ] = Column( mag_type2.T )   
t['REDDENING'         ] = Column( data[1].data['ebv_galaxy'][eRO_DE])
t['DATE_EARLIEST'     ] = Column( n.ones(N_tot)*59215. )  
t['DATE_LATEST'       ] = Column( n.ones(N_tot)*66520 )   

print( path_2_out )
if os.path.isfile(path_2_out):
	os.system("rm "+path_2_out)

t.write(path_2_out, format='fits')


s3 = ( data_3[1].data['DEC']>-80)&( data_3[1].data['DEC']<5)&(data_3[1].data['R_MAG']<23.5) # 23.2
N_tot = len(data_3[1].data['RA'][s3])

s3_name = 'AGN_IR'
sub_cats = n.zeros(N_tot).astype('U20')
sub_cats[:] = s3_name

ruleset_array = n.zeros(N_tot).astype('U20')
ruleset_array[:] = "AGN_ALL_3PC"

r_mag = data_3[1].data['R_MAG'][s3]

#AGN_ALL_*PC_rXXY = ruleset to be used for all targets with r.gt.22.0, with magnitude XX.Y<r<=XX.Y+0.1
#faint_ruleset = n.array(['AGN_ALL_3PC_r'+str(int(rmag*10)) for rmag in r_mag ])
#ruleset_array[r_mag>22] = faint_ruleset[r_mag>22]

mag_type = n.zeros(N_tot).astype('U9')
mag_type[:] = "SDSS_r_AB"
mag_type2 = n.zeros((1,N_tot)).astype('U9')
mag_type2[0]=mag_type


list_2_check = n.unique(data_3[1].data['TEMPLATE'][s3])
for el in list_2_check:
	print(os.path.isfile(os.path.join('templates',el)), el)

t = Table()

t['NAME'              ] = Column( n.arange(N_tot).astype('str') )
t['RA'                ] = Column( data_3[1].data['RA'] [s3] )   
t['DEC'               ] = Column( data_3[1].data['DEC']    [s3] )   
t['PMRA'              ] = Column( n.zeros(N_tot) )   
t['PMDEC'             ] = Column( n.zeros(N_tot) )   
t['EPOCH'             ] = Column( n.ones(N_tot)*2000.0 )   
t['RESOLUTION'        ] = Column( n.ones((1,N_tot)).T ) 
t['SUBSURVEY'         ] = Column( sub_cats )  
t['PRIORITY'          ] = Column( n.ones(N_tot)*100.0 )
t['TEMPLATE'          ] = Column( data_3[1].data['TEMPLATE'][s3] )
t['RULESET'           ] = Column( ruleset_array ) 
t['REDSHIFT_ESTIMATE' ] = Column( data_3[1].data['REDSHIFT'][s3].astype('float') )
t['REDSHIFT_ERROR'    ] = Column( n.ones(N_tot)*0.00001 )
t['EXTENT_PARAMETER'  ] = Column( n.ones(N_tot) )   
t['MAG'               ] = Column( r_mag )   
t['MAG_ERR'           ] = Column( n.ones(N_tot)*0.01 )   
t['MAG_TYPE'          ] = Column( mag_type2.T )   
t['REDDENING'         ] = Column( data_3[1].data['EBV'][s3] )
t['DATE_EARLIEST'     ] = Column( n.ones(N_tot)*59215. )  
t['DATE_LATEST'       ] = Column( n.ones(N_tot)*66520 )   


print( path_2_out_3 )
if os.path.isfile(path_2_out_3):
	os.system("rm "+path_2_out_3)

t.write(path_2_out_3, format='fits')

# now concatenates the tables
t1=fits.open(path_2_out)[1].data
t2=fits.open(path_2_out_3)[1].data
N_tot = len(t1)+len(t2)

mag_type = n.zeros(N_tot).astype('U9')
mag_type[:] = "SDSS_r_AB"
mag_type2 = n.zeros((1,N_tot)).astype('U9')
mag_type2[0][:]="SDSS_r_AB"

t_o = Table()
t_o['NAME'] = Column( n.arange(N_tot).astype('str') )
for cname in t1.columns.names[1:]:
	print(cname)
	if cname=='MAG_TYPE':
		print('mag type')
		t_o[cname] = Column( mag_type2.T, dtype='U9' ) 
	else:
		t_o[cname] = Column( n.hstack((t1[cname], t2[cname])) ) 

path_final_file="S6_summary.fits"
if os.path.isfile(path_final_file):
	os.system("rm "+path_final_file)
t_o.write(path_final_file, format='fits')

#stilts_cmd = 'java -jar ~/software/stilts.jar tcat ifmt=fits in='+path_2_out_3+' in='+path_2_out+' omode=out out=S6_summary_may2nd.fits'
#os.system(stilts_cmd)

test = fits.open(path_final_file)

test2 = fits.open(path_2_out_3)