"""
Creates a fits catalog containing the 4FS input columns.

Create the ID array
remove unwanted column
check rulesets and templates names are correct

005_4_add_footprint.py


python 005_4_add_footprint.py S5_BCG_4MOST.fit
python 005_4_add_footprint.py LyA_4MOST.fits
python 005_4_add_footprint.py BG_4MOST.fits
python 005_4_add_footprint.py LRG_4MOST.fits
python 005_4_add_footprint.py QSO_4MOST.fits
python 005_4_add_footprint.py AGN_WIDE_4MOST.fits
python 005_4_add_footprint.py AGN_DEEP_4MOST.fits
python 005_4_add_footprint.py FILAMENTS5_4MOST.fits
python 005_4_add_footprint.py S5_CGAL_4MOST.fit
python 005_4_add_footprint.py AGN_IR_4MOST.fits
python 005_4_add_footprint.py ELG_4MOST.fits

cd $MD10
gzip -k --rsyncable  $MD10/S5_BCG_4MOST.fit
gzip -k --rsyncable  $MD10/LyA_4MOST.fits
gzip -k --rsyncable  $MD10/BG_4MOST.fits
gzip -k --rsyncable  $MD10/S5_CGAL_4MOST.fit
gzip -k --rsyncable  $MD10/AGN_IR_4MOST.fits
gzip -k --rsyncable  $MD10/LRG_4MOST.fits
gzip -k --rsyncable  $MD10/QSO_4MOST.fits
gzip -k --rsyncable  $MD10/AGN_WIDE_4MOST.fits
gzip -k --rsyncable  $MD10/AGN_DEEP_4MOST.fits
gzip -k --rsyncable  $MD10/FILAMENTS5_4MOST.fits

rsync $MD10/S5_BCG_4MOST.fit.gz      $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/LyA_4MOST.fits.gz        $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/BG_4MOST.fits.gz         $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/LRG_4MOST.fits.gz        $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/S5_CGAL_4MOST.fit.gz     $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/AGN_IR_4MOST.fits.gz     $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/QSO_4MOST.fits.gz        $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/AGN_WIDE_4MOST.fits.gz   $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/AGN_DEEP_4MOST.fits.gz   $HOME/wwwDir/MultiDark/MD10/
rsync $MD10/FILAMENTS5_4MOST.fits.gz $HOME/wwwDir/MultiDark/MD10/

gzip -k --rsyncable  $MD10/ELG_4MOST.fits
rsync $MD10/ELG_4MOST.fits.gz        $HOME/wwwDir/MultiDark/MD10/


topcat -stilts plot2plane \
   xpix=1376 ypix=586 \
   xlabel=REDSHIFT_ESTIMATE ylabel=galaxy_stellar_mass fontsize=18 \
   xmin=0.001 xmax=3 ymin=8 ymax=12.53 \
   legend=true legpos=1.0,0.0 \
   nlevel=6 smooth=18 scaling=log \
   layer_01=Contour \
      in_01=/data/data/MultiDark/MD_1.0Gpc/AGN_WIDE_4MOST.fits \
      x_01=REDSHIFT_ESTIMATE y_01=galaxy_stellar_mass \
      zero_01=1.007 \
      leglabel_01='AGN eROSITA' \
   layer_02=Contour \
      in_02=/data/data/MultiDark/MD_1.0Gpc/QSO_4MOST.fits \
      x_02=REDSHIFT_ESTIMATE y_02=galaxy_stellar_mass \
      color_02=orange \
      leglabel_02='QSO S8' \
   layer_03=Contour \
      in_03=/data/data/MultiDark/MD_1.0Gpc/AGN_IR_4MOST.fits \
      x_03=REDSHIFT_ESTIMATE y_03=galaxy_stellar_mass \
      color_03=cyan \
      leglabel_03='AGN IR' \
   layer_04=Contour \
      in_04=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/S5GAL_000004.fit \
      x_04=Z y_04=Mstar \
      color_04=magenta \
      leglabel_04='S5 filament' \
   layer_05=Contour \
      in_05=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/BG_000004.fit \
      x_05=Z y_05=Mstar \
      color_05=blue \
      leglabel_05='S8 BG' \
   layer_06=Contour \
      in_06=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/LRG_000004.fit \
      x_06=Z y_06=Mstar \
      color_06=green \
      leglabel_06='S8 LRG' \
   layer_07=Contour \
      in_07=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/ELG_000004.fit \
      x_07=Z y_07=Mstar \
      color_07=grey \
      leglabel_07='S8 ELG' \
   layer_08=Contour \
      in_08=/data/data/MultiDark/MD_1.0Gpc/LyA_4MOST.fits \
      x_08=REDSHIFT_ESTIMATE y_08=galaxy_stellar_mass \
      color_08=black \
      leglabel_08='S8 Lya' \
   layer_09=Contour \
      in_09=/data/data/MultiDark/MD_1.0Gpc/S5_CGAL_4MOST.fit \
      x_09=redshift_R y_09=galaxy_stellar_mass \
      color_09=light_grey \
      leglabel_09='S5 cluster GAL' \
   layer_10=Contour \
      in_10=/data/data/MultiDark/MD_1.0Gpc/S5_BCG_4MOST.fit \
      x_10=redshift_R y_10=galaxy_stellar_mass \
      color_10=yellow \
      leglabel_10='S5 BCG' \
    omode=out out=/home/comparat/data/MultiDark/MD_1.0Gpc/stellarMass-redshift.png


topcat -stilts plot2plane \
   xpix=1376 ypix=586 \
   xlabel=REDSHIFT_ESTIMATE ylabel='MAG / mag' fontsize=18 \
   xmin=0.001 xmax=3 ymin=14 ymax=25 \
   legend=true legpos=1.0,0.0 \
   nlevel=6 smooth=18 scaling=log \
   layer_01=Contour \
      in_01=/data/data/MultiDark/MD_1.0Gpc/AGN_WIDE_4MOST.fits \
      x_01=REDSHIFT_ESTIMATE y_01=MAG \
      zero_01=1.007 \
      leglabel_01='AGN eROSITA' \
   layer_02=Contour \
      in_02=/data/data/MultiDark/MD_1.0Gpc/QSO_4MOST.fits \
      x_02=REDSHIFT_ESTIMATE y_02=MAG \
      color_02=orange \
      leglabel_02='QSO S8' \
   layer_03=Contour \
      in_03=/data/data/MultiDark/MD_1.0Gpc/AGN_IR_4MOST.fits \
      x_03=REDSHIFT_ESTIMATE y_03=MAG \
      color_03=cyan \
      leglabel_03='AGN IR' \
   layer_04=Contour \
      in_04=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/S5GAL_000004.fit \
      x_04=Z y_04=rfib \
      color_04=magenta \
      leglabel_04='S5 filament' \
   layer_05=Contour \
      in_05=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/BG_000004.fit \
      x_05=Z y_05=rfib \
      color_05=blue \
      leglabel_05='S8 BG' \
   layer_06=Contour \
      in_06=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/LRG_000004.fit \
      x_06=Z y_06=rfib \
      color_06=green \
      leglabel_06='S8 LRG' \
   layer_07=Contour \
      in_07=/data/data/MultiDark/MD_1.0Gpc/cat_SHAM_COSMO/ELG_000004.fit \
      x_07=Z y_07=rfib \
      color_07=grey \
      leglabel_07='S8 ELG' \
   layer_08=Contour \
      in_08=/data/data/MultiDark/MD_1.0Gpc/LyA_4MOST.fits \
      x_08=REDSHIFT_ESTIMATE y_08=MAG \
      color_08=black \
      leglabel_08='S8 Lya' \
   layer_09=Contour \
      in_09=/data/data/MultiDark/MD_1.0Gpc/S5_CGAL_4MOST.fit \
      x_09=redshift_R y_09=MAG \
      color_09=light_grey \
      leglabel_09='S5 cluster GAL' \
   layer_10=Contour \
      in_10=/data/data/MultiDark/MD_1.0Gpc/S5_BCG_4MOST.fit \
      x_10=redshift_R y_10=MAG \
      color_10=yellow \
      leglabel_10='S5 BCG' \
    omode=out out=/home/comparat/data/MultiDark/MD_1.0Gpc/fibermagnitude-redshift.png

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
from matplotlib.path import Path 
print('adds survey bits')
print('------------------------------------------------')
print('------------------------------------------------')

MASKDIR = os.path.join(os.environ['GIT_AGN_MOCK'],'data', 'masks')

if len(sys.argv)==1:
	dir_2_MOCK = sys.argv[1]
if len(sys.argv)==2:
	env = sys.argv[1] 
	file_name = sys.argv[2]
	root_dir = os.path.join(os.environ[env])
	dir_2_MOCK  = os.path.join(root_dir, file_name)

nl = lambda selection : len(selection.nonzero()[0])


import numpy as np
from astropy import units
from astropy.coordinates import SkyCoord
import pymangle


# get survey footprint bitlist
def fourmost_get_survbitlist():
	mydict                = {}
	mydict['des']         = 0
	mydict['kidss']       = 1
	mydict['kidsn']       = 2
	mydict['atlassgcnotdes'] = 3
	mydict['atlasngc']    = 4
	mydict['kabs']        = 5
	mydict['vhsb10']      = 6
	mydict['vhsb15']      = 7
	mydict['vhsb20']      = 8
	mydict['vhsb20clean'] = 9
	mydict['desi']        = 10
	mydict['erosita']     = 11
	mydict['waveswide']   = 12
	mydict['euclid']      = 13
	mydict['s8elg']       = 14
	mydict['s8']          = 15
	return mydict


# 4most s8 footprint
def fourmost_get_s8foot(ra,dec):
	svbdict = fourmost_get_survbitlist()
	survbit = np.zeros(len(ra),dtype=int)
	for bname in ['des','desi','kidss','kidsn','atlassgcnotdes','atlasngc']:
		isb           = fourmost_get_survbit_indiv(ra,dec,bname)
		survbit[isb] += 2**svbdict[bname]
	# bg/lrgi/qso/lya
	iss8  = (
		(( (survbit & 2**svbdict['des'])>0)   & ( (survbit & 2**svbdict['desi'])==0)) | # des not desi
		(( (survbit & 2**svbdict['kidss'])>0) | (((survbit & 2**svbdict['kidsn'])>0) & (ra>154))) | # kids
		((((survbit & 2**svbdict['atlassgcnotdes'])>0) | ( (survbit & 2**svbdict['atlasngc'])>0)) & ((survbit & 2**svbdict['desi'])==0) & (dec<-10)) # atlas not desi
	)
	# elg
	iss8elg = (((ra>330) | (ra<90)) & (dec>-35.5) & (dec<-26)) | ((ra>50) & (ra<90) & (dec>-40) & (dec<-26))
	return iss8,iss8elg



# get survey footprint bit
def fourmost_get_survbit_indiv(ra,dec,bname):
	#
	# Galactic l,b
	c       = SkyCoord(ra=ra, dec=dec, frame='fk5')
	l,b     = c.galactic.l.value,c.galactic.b.value
	lon,lat = c.barycentrictrueecliptic.lon.degree,c.barycentrictrueecliptic.lat.degree
	# 4most/s8 , s8elg
	if (bname in ['s8elg','s8']):
		iss8,iss8elg = fourmost_get_s8foot(ra,dec)
		if (bname=='s8elg'):
			keep = iss8elg
		else: # s8
			keep = iss8
	# desi: no ply file for the moment...
	elif (bname=='desi'):
		# desi sgc
		polyra = np.array([0 ,-25,-35,-50,-54,-45,10,10,  60, 70, 70,53,42,42,38,0])
		polydec= np.array([33,33, 25,  8,  -8, -15,-15,-20,-20,-15,0, 0, 10,20,33,33])
		sgcpoly = Path(np.concatenate(
							(polyra. reshape((len(polyra),1)),
							 polydec.reshape((len(polyra),1))),
							axis=1))
		# desi ngc
		polyra = np.array([275,107,115,130,230,230,230,255,265,275])
		polydec= np.array([33, 33, 12, -10,-10, -2, -2, -2,  13, 33])
		ngcpoly = Path(np.concatenate(
							(polyra. reshape((len(polyra),1)),
							 polydec.reshape((len(polyra),1))),
							axis=1))
		#
		tmpradec         = np.transpose(np.array([ra,dec]))
		tmp              = (ra>300)
		tmpradec[tmp,0] -= 360.
		keep = np.zeros(len(ra),dtype=bool)
		for poly in [sgcpoly,ngcpoly]:
			keep[poly.contains_points(tmpradec)] = True
	elif (bname=='erosita'):
		# johan email 30/07/2018 14:12
		keep = (abs(b)>15) & (l>180)
	elif (bname=='waveswide'):
		# https://wavesurvey.org/project/survey-design/
		keep = (((ra>155) & (ra<240) & (dec>-5) & (dec<5))
				|
				(((ra>330) | (ra<50)) & (dec>-36) & (dec<-26)))
	elif (bname=='euclid'):
		keep = (np.abs(b)>=30) & (np.abs(lat)>5.)
	else:
		if (bname[:3]=='vhs'):
			## -70<dec<0
			mng    = pymangle.Mangle(MASKDIR+'/vhsdec.ply')
			polyid = mng.polyid(ra,dec)
			keepdec= (polyid!=-1)
			# |b|>bmin
			mng    = pymangle.Mangle(MASKDIR+'/'+bname[3:6]+'.ply')
			polyid = mng.polyid(l,b)
			keepb  = (polyid!=-1)
			##
			keep   = (keepdec) & (keepb)
		else:
			mng    = pymangle.Mangle(MASKDIR+'/'+bname+'.ply')
			polyid = mng.polyid(ra,dec)
			keep   = (polyid!=-1)
		if (bname=='vhsb20clean'):
			## Jext<0.1 and low nstar selection [both VHS and DES]
			ra60    = (ra>55)  & (ra<65)  & (dec>-5)  & (dec<0)
			ra70    = (ra>67)  & (ra<72)  & (dec>-16) & (dec<-13)
			ra100   = (ra>50)  & (ra<180) & (dec>-20) & (dec<0)   & (b>-23) & (b<0)
			ra120   = (ra>100) & (ra<180)                         & (b>-15) & (b<0)
			ra230   = (ra>228) & (ra<270) & (dec>-40) & (dec<-20) & (b>0) 
			ra250   = (ra>235) & (ra<270) & (dec>-20) & (dec<0)   & (b>0)
			ra300   = (ra>180) & (ra<360) & (dec>-70) & (dec<0)   & (b>-25) & (b<0)
			LMC     = (ra>70)  & (ra<90)  & (dec>-70) & (dec<-65)
			keep    = ((keep) & 
						(~ra60)  & (~ra70)  & 
						(~ra100) & (~ra120) & 
						(~ra230) & (~ra250) & (~ra300) &
						(~LMC))
	print( bname, len(ra[keep]))
	return keep


def fourmost_get_survbit(ra,dec):
	bdict   = fourmost_get_survbitlist()
	survbit = np.zeros(len(ra),dtype='int')
	for bname in bdict.keys():
		print(bname)
		isb           = fourmost_get_survbit_indiv(ra,dec,bname)
		survbit[isb] += 2**bdict[bname]
	return survbit


t_survey  = Table.read(dir_2_MOCK  )
print(dir_2_MOCK)
mask_bit = fourmost_get_survbit(t_survey['RA'],t_survey['DEC'])
print(t_survey.columns.keys())
if 'MASK_BIT' in t_survey.columns.keys():
	t_survey['MASK_BIT'] = mask_bit
else:
	t_survey.add_column(Column(name='MASK_BIT'  ,data=mask_bit, unit=''))
	
t_survey.write (dir_2_MOCK  , overwrite=True)

