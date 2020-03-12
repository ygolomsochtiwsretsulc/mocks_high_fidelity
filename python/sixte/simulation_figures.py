"""
##!bin/bash

#topcat -stilts plot2plane xpix=400 ypix=400 layer=Mark in=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100.fits icmd=explodeall icmd='select " SRC_ID_1>=1e9"' fontsize=16 fontweight=bold xmin=206 xmax=218 ymin=-54 ymax=-42 xlabel='RA / deg' ylabel='DEC / deg' x=RA y=DEC shading=transparent size=0 color=black opaque=20 omode=out out=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100_A.png

#topcat -stilts plot2plane xpix=400 ypix=400 layer=Mark in=/home/comparat/data/erosim/eRASS/pointSourceCat-v0.1/rawevt_100.fits icmd=explodeall icmd='select " SRC_ID_1>=1e9"' fontsize=16 fontweight=bold xmin=206 xmax=218 ymin=-54 ymax=-42 xlabel='RA / deg' ylabel='DEC / deg' x=RA y=DEC shading=transparent size=0 color=black opaque=10 omode=out out=/home/comparat/data/erosim/eRASS/pointSourceCat-v0.1/rawevt_100_A.png



#topcat -stilts plot2plane xpix=400 ypix=400 layer=Mark in=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100.fits icmd=explodeall icmd='select "SRC_ID_1>=1e6 && SRC_ID_1<1e9"' fontsize=16 fontweight=bold xmin=206 xmax=218 ymin=-54 ymax=-42 xlabel='RA / deg' ylabel='DEC / deg' x=RA y=DEC shading=transparent size=0 color=black opaque=10 omode=out out=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100_C.png

#topcat -stilts plot2plane xpix=400 ypix=400 layer=Mark in=/home/comparat/data/erosim/eRASS/pointSourceCat-v0.1/rawevt_100.fits icmd=explodeall icmd='select "SRC_ID_1>=1e6 && SRC_ID_1<1e9"' fontsize=16 fontweight=bold xmin=206 xmax=218 ymin=-54 ymax=-42 xlabel='RA / deg' ylabel='DEC / deg' x=RA y=DEC shading=transparent size=0 color=black opaque=10 omode=out out=/home/comparat/data/erosim/eRASS/pointSourceCat-v0.1/rawevt_100_C.png


#topcat -stilts plot2plane \
   #xpix=1200 ypix=600 fontsize=16 fontweight=bold \
   #xlabel='RA / deg' ylabel='DEC / deg' \
   #xmin=210 xmax=214 ymin=-50 ymax=-48 \
   #legend=true legpos=1.0,1.0 \
   #x=RA y=DEC \
   #layer_1=Mark \
      #in_1=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100.fits \
      #shading_1=transparent size_1=0 color_1=grey opaque_1=8 \
   #layer_2=Mark \
      #in_2=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100.fits \
       #icmd_2=explodeall \
       #icmd_2='select SRC_ID_1>0' \
      #shading_2=transparent size_2=0 color_2=grey opaque_2=6 \
   #layer_3=Mark \
      #in_3=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/rawevt_100.fits \
       #icmd_3=explodeall \
       #icmd_3='select SRC_ID_1>1e6' \
      #shading_3=transparent size_3=0 color_3=black opaque_3=4 \
      #leglabel_3='eRASS8 events' \
   #layer_4=Mark \
      #in_4=/home/comparat/data/erosim/eRASS/eSASS-v11-erass8-jc/020_MLCat.fits \
      #shading_4=auto shape_4=open_square size_4=5 color_4=blue \
      #leglabel_4='eRASS8 catalogue' \
   #layer_5=Mark \
      #in_5=/home/comparat/data/erosim/eRASS/pointSourceCat-v0.1/020_MLCat.fits \
      #shading_5=auto shape_5=open_circle size_5=5 \
      #leglabel_5='eRASS2 catalogue' \
   #legseq=_3,_4,_5 \
   #omode=out out=/home/comparat/data/erosim/eRASS/detection.png

#Preparation for the plotting
#!/bin/bash 

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd1_evt.fits > list_ccd1.lis
stilts tcat in=@list_ccd1.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd1.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd2_evt.fits > list_ccd2.lis
stilts tcat in=@list_ccd2.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd2.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd3_evt.fits > list_ccd3.lis
stilts tcat in=@list_ccd3.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd3.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd4_evt.fits > list_ccd4.lis
stilts tcat in=@list_ccd4.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd4.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd5_evt.fits > list_ccd5.lis
stilts tcat in=@list_ccd5.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd5.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd6_evt.fits > list_ccd6.lis
stilts tcat in=@list_ccd6.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd6.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/???/erass_ccd7_evt.fits > list_ccd7.lis
stilts tcat in=@list_ccd7.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd7.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccd* > list_ccdA.lis
stilts tcat in=@list_ccdA.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits ofmt=fits

mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_1
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_2
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_3
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_4
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_5
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_6
mkdir /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cam_7

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.003 auxmax=0.3 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cluster_tAll.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.003 auxmax=0.3 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000/2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS1' \
      omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cluster_eRASS1.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.003 auxmax=0.3 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS2' \
      omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cluster_eRASS2.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.003 auxmax=0.3 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000*2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS4' \
      omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cluster_eRASS4.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.003 auxmax=0.3 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000*3' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS6' \
      omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cluster_eRASS6.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.003 auxmax=0.3 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000*4' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS8' \
      omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/cluster_eRASS8.png

cp -r /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures $GIT_AGN_MOCK/figures/eRASS8_cluster_MD40/

wwwDir/erosita_stuff/
#cd data/eRoMok
#wget http://www.mpe.mpg.de/~comparat/erosita_stuff/simulated_photons.fits

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 1

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 2

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 3

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 4

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 5

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 6

ds52
cd lss_mock_dev/python/sixte/
pyCONDA
python simulation_figures.py 7


cd /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures/

convert -delay 20 -loop 0 cam_1/agn_*.png cam_1/agn_year1.gif
convert -delay 20 -loop 0 cam_2/agn_*.png cam_2/agn_year1.gif
#convert -delay 20 -loop 0 cam_3/agn_*.png cam_3/agn_year1.gif
#convert -delay 20 -loop 0 cam_4/agn_*.png cam_4/agn_year1.gif
#convert -delay 20 -loop 0 cam_5/agn_*.png cam_5/agn_year1.gif
#convert -delay 20 -loop 0 cam_6/agn_*.png cam_6/agn_year1.gif
convert -delay 20 -loop 0 cam_7/agn_*.png cam_7/agn_year1.gif

convert -delay 20 -loop 0 cam_1/cluster_*.png cam_1/cluster_year1.gif
convert -delay 20 -loop 0 cam_2/cluster_*.png cam_2/cluster_year1.gif
#convert -delay 20 -loop 0 cam_3/cluster_*.png cam_3/cluster_year1.gif
#convert -delay 20 -loop 0 cam_4/cluster_*.png cam_4/cluster_year1.gif
#convert -delay 20 -loop 0 cam_5/cluster_*.png cam_5/cluster_year1.gif
#convert -delay 20 -loop 0 cam_6/cluster_*.png cam_6/cluster_year1.gif
convert -delay 20 -loop 0 cam_7/cluster_*.png cam_7/cluster_year1.gif

"""

import numpy as n
import os, sys
from astropy.table import Table
N_cam = sys.argv[1] # '1'
root_dir = '/data40s/erosim/eRASS/eRASS8'
path_2_file = os.path.join(root_dir, "simulated_photons_ccd"+N_cam+".fits")
fig_dir = os.path.join(root_dir, 'figures', 'cam_' + N_cam)

data = Table.read(path_2_file)
agn = (data['SRC_ID_1']>1e9)
cluster = (data['SRC_ID_1']<1e9)
#cluster_ids = data['SRC_ID_1'][]
#n.unique()

def plot_sim_cluster(time_max):
	str_t = str(int(time_max)).zfill(8)
	str_t_label = str(int(time_max/(3600*24.))).zfill(4)
	c1="topcat -stilts plot2sky xpix=1024 ypix=512 projection=aitoff legend=true legpos=0.9,0.9 layer=Mark "
	c2="""in="""+path_2_file+""" icmd='select "TIME<"""+str_t+""" && SRC_ID_1<1e9"' """
	c3="""lon=RA lat=DEC shading=transparent size=0 color=black opaque=200 leglabel="""+str_t_label
	c4="""days omode=out out="""+fig_dir+"""/cluster_t"""+str_t+""".png"""
	return c1+c2+c3+c4

for tt in n.arange(0.25,41,0.25):
	cmd = plot_sim_cluster(tt*1e6)
	print(cmd)
	os.system(cmd)

def plot_sim_agn(time_max):
	str_t = str(int(time_max)).zfill(8)
	str_t_label = str(int(time_max/(3600*24.))).zfill(4)
	c1="topcat -stilts plot2sky xpix=1024 ypix=512 projection=aitoff legend=true legpos=0.9,0.9 layer=Mark "
	c2="""in="""+path_2_file+""" icmd='select "TIME<"""+str_t+""" && SRC_ID_1>1e9"' """
	c3="""lon=RA lat=DEC shading=transparent size=0 color=black opaque=500 leglabel="""+str_t_label
	c4="""days omode=out out="""+fig_dir+"""/agn_t"""+str_t+""".png"""
	return c1+c2+c3+c4

for tt in n.arange(0.25,41,0.25):
	cmd = plot_sim_agn(tt*1e6)
	print(cmd)
	os.system(cmd)


