#!bin/bash

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd1_evt.fits > list_ccd1.lis
stilts tcat in=@list_ccd1.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd1.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd2_evt.fits > list_ccd2.lis
stilts tcat in=@list_ccd2.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd2.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd3_evt.fits > list_ccd3.lis
stilts tcat in=@list_ccd3.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd3.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd4_evt.fits > list_ccd4.lis
stilts tcat in=@list_ccd4.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd4.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd5_evt.fits > list_ccd5.lis
stilts tcat in=@list_ccd5.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd5.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd6_evt.fits > list_ccd6.lis
stilts tcat in=@list_ccd6.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd6.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd7_evt.fits > list_ccd7.lis
stilts tcat in=@list_ccd7.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd7.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd?.fits > list_ccdA.lis
stilts tcat in=@list_ccdA.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits ofmt=fits


###########################################################
###########################################################
# CLUSTER
###########################################################
###########################################################

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      omode=out out=$GIT_AGN_MOCK/figures/MD40/sixte/cluster_tAll.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000/2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS1' \
      omode=out out=$GIT_AGN_MOCK/figures/MD40/sixte/cluster_eRASS1.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS2' \
      omode=out out=$GIT_AGN_MOCK/figures/MD40/sixte/cluster_eRASS2.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000*2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS4' \
      omode=out out=$GIT_AGN_MOCK/figures/MD40/sixte/cluster_eRASS4.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000*3' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS6' \
      omode=out out=$GIT_AGN_MOCK/figures/MD40/sixte/cluster_eRASS6.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits icmd='select TIME<31536000*4' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS8' \
      omode=out out=$GIT_AGN_MOCK/figures/MD40/sixte/cluster_eRASS8.png

###########################################################
###########################################################
# CLUSTER
###########################################################
###########################################################

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/cluster_tAll.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000/2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS1' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/cluster_eRASS1.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS2' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/cluster_eRASS2.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000*2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS4' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/cluster_eRASS4.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000*3' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS6' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/cluster_eRASS6.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000*4' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS8' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/cluster_eRASS8.png

###########################################################
###########################################################
# AGN
###########################################################
###########################################################

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccdA.fits \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/agn_tAll.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000/2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS1' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/agn_eRASS1.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS2' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/agn_eRASS2.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000*2' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS4' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/agn_eRASS4.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000*3' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS6' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/agn_eRASS6.png

topcat -stilts plot2sky \
   xpix=1543 ypix=721 \
   projection=aitoff crowd=7.2 gridaa=true texttype=antialias fontsize=16 fontweight=bold \
   auxmap=accent auxfunc=log auxmin=0.01 auxmax=10. auxcrowd=3.5 \
   auxvisible=true auxlabel='density of counts per arcmin**2' \
   legend=false \
   layer=SkyDensity \
      in=$HOME/erosim/eRASS/eRASS8_agn_MD10/simulated_photons_ccdA.fits icmd='select TIME<31536000*4' \
      lon=RA lat=DEC \
      level=5 combine=count-per-unit perunit=arcmin2 \
      leglabel='1: eRASS8' \
      omode=out out=$GIT_AGN_MOCK/figures/MD10/sixte/agn_eRASS8.png


      

wwwDir/erosita_stuff/
cd data/eRoMok
wget http://www.mpe.mpg.de/~comparat/erosita_stuff/simulated_photons.fits

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


cd $HOME/erosim/eRASS/eRASS8_cluster_MD10/figures/

convert -delay 20 -loop 0 cam_1/agn_*.png cam_1/agn_year1.gif
convert -delay 20 -loop 0 cam_2/agn_*.png cam_2/agn_year1.gif
convert -delay 20 -loop 0 cam_3/agn_*.png cam_3/agn_year1.gif
convert -delay 20 -loop 0 cam_4/agn_*.png cam_4/agn_year1.gif
convert -delay 20 -loop 0 cam_5/agn_*.png cam_5/agn_year1.gif
convert -delay 20 -loop 0 cam_6/agn_*.png cam_6/agn_year1.gif
convert -delay 20 -loop 0 cam_7/agn_*.png cam_7/agn_year1.gif

convert -delay 20 -loop 0 cam_1/cluster_*.png cam_1/cluster_year1.gif
convert -delay 20 -loop 0 cam_2/cluster_*.png cam_2/cluster_year1.gif
convert -delay 20 -loop 0 cam_3/cluster_*.png cam_3/cluster_year1.gif
convert -delay 20 -loop 0 cam_4/cluster_*.png cam_4/cluster_year1.gif
convert -delay 20 -loop 0 cam_5/cluster_*.png cam_5/cluster_year1.gif
convert -delay 20 -loop 0 cam_6/cluster_*.png cam_6/cluster_year1.gif
convert -delay 20 -loop 0 cam_7/cluster_*.png cam_7/cluster_year1.gif

