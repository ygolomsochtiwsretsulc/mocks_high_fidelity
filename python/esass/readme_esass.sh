 ... offending file: /home/erosita/caldb/data/erosita/tm1/caldb.indx
 erbox: eSASSusers_200302 Mar 13 17:44:18 2020
 erbox/ero_cifsl: accessing cal file tm1_2dpsf_190220v03.fits
 Filtered source list:        36556
 erbox/ero_cifsl: **ERROR3** Unable to open CIF
 erbox/ero_cifsl: **STOP** Fatal - aborting
 erbox: FAILED
 erbackmap: eSASSusers_200302 Mar 13 17:44:18 2020
 erbackmap/FGET_BOXLIST: **ERROR3** error opening file /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/016/016_
 erbackmap: FAILED
 erbox: eSASSusers_200302 Mar 13 17:44:18 2020
 erbox/fget_r4compress: **ERROR3** Error opening FITS file
 erbox/ERBOX_IN: **ERROR3** Cannot read background map
 erbox: FAILED
 erbackmap: eSASSusers_200302 Mar 13 17:44:18 2020
 erbackmap/FGET_BOXLIST: **ERROR3** error opening file /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/016/016_
 erbackmap: FAILED
outprefix 016_

Traceback (most recent call last):
  File "/home/comparat/software/linux/lss_mock_dev/python/esass/eRO_pipe_det.py", line 292, in <module>
    assert all(os.path.isfile(bm) for bm in BkgMapFiles[Nmlin])
AssertionError


ds52: no
ds43

radec2xy file=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/104/tmp_104.fits ra0=101.24999999999999 dec0=35.68533471265204

chgrp erosim 0??/*
chgrp erosim 1??/*
chgrp erosim 2??/*
chgrp erosim 3??/*
chgrp erosim 4??/*
chgrp erosim 5??/*
chgrp erosim 6??/*
chgrp erosim 7??/*

source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh

cd $GIT_AGN_MOCK/python/esass
ls /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/???_0216_Sc1Cat.fits > Sc1Cat_done.list 


# rm /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/*0216*
# rm /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/*.log

First test (Georg) : 
1. eRO_pipe_images

- single band: 0.2, 2.2
- single band: 0.3, 2.2
- single band: 0.2, 2.3
- single band: 0.3, 2.3
- single band: 0.5, 2.3 
- single image 0.2-0.5 keV + 0.6-2.3 keV to exclude a feature at 0.5-0.6 keV (strong background feature to be singled out)
  + find accurately the position and width of the background peak to define the 2 bands

2. eRO_pipe_det 

Change the last erbox and erbackmap parameters, not the first 

2.1. erbox, line 155, impacts background calculation
 - likemin = 5, 6, 7
 - nruns = 2, 3

2.2. er_backmap
 - fitmethod=smooth
 - smoothval = 15.  
 - snr = 10, 20, 30, 40, 50, 60 (increase, then BKG becomes smoother) Test the values
 - excesssigma = 1000.  
 - nfitrun = 3 

2.3. ermldet, line 106, parameters 
 - extentmodel = beta
 - ERML_CUTRAD    = 12, 15, 20 in pixels (4 arc seconds per pixels) 
 - ERML_MULTRAD   = 12, 15, 20    (where it starts to fit multiple sources)
 - ERML_likemin   =  5
 - ERML_extlikemin=  5
 - ERML_extmin    =  2
 - ERML_extmax    = 12, 15, 20 not higher than ERML_CUTRAD

2.4
 - Do not execute the second ermldet. Only there to improve point source detection
 - Do not execute ApeTool
# are inputs there: yes !
grep 000 $GIT_AGN_MOCK/python/*/remaining*.sh
grep 699 $GIT_AGN_MOCK/python/*/remaining*.sh
grep 767 $GIT_AGN_MOCK/python/*/remaining*.sh
grep 709 $GIT_AGN_MOCK/python/*/remaining*.sh
grep 608 $GIT_AGN_MOCK/python/*/remaining*.sh
grep 528 $GIT_AGN_MOCK/python/*/remaining*.sh

grep 000 $GIT_AGN_MOCK/python/*/esass_commands_detection.sh
grep 699 $GIT_AGN_MOCK/python/*/esass_commands_detection.sh
grep 767 $GIT_AGN_MOCK/python/*/esass_commands_detection.sh
grep 709 $GIT_AGN_MOCK/python/*/esass_commands_detection.sh
grep 608 $GIT_AGN_MOCK/python/*/esass_commands_detection.sh
grep 528 $GIT_AGN_MOCK/python/*/esass_commands_detection.sh


nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000 -V 13 > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_000_eSASSusers_200302.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699 -V 13 > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_699_eSASSusers_200302.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767 -V 13 > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_767_eSASSusers_200302.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709 -V 13 > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_709_eSASSusers_200302.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608 -V 13 > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_608_eSASSusers_200302.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528 -V 13 > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_528_eSASSusers_200302.log & 


cd $GIT_AGN_MOCK/python/esass
# generate this file 
# healpix_radius_nested_all.dat
# with 
python tabulate_center.py

# write the commands to convert all event files into a single one to be fed into eSASS
python calibrate_events_WriteCommands.py
# it creates 
# a set of 768 calibrate_events.py commands written in 
calibrate_events.sh
# run all the commands inside this script after having sourced eSASS
# it needs evtool eSASS command (i.e. to source eSASS)

python calibrate_events_WriteCommands_ready_to_go.py

nohup sh calibrate_events_ready_to_go.sh > logs/calibrate_events_list_ready_4_go.log & 
nohup sh calibrate_events_ready_to_go.sh > logs/calibrate_events_list_ready_5_go.log & 

nohup sh calibrate_events_list_by_hand.sh > logs/calibrate_events_list_by_hand.log & 

nohup sh calibrate_events_list_by_hand.sh > logs/calibrate_events_list_by_hand.log & 
nohup sh calibrate_events.sh > logs/calibrate_events_all.log & 

cd /home/comparat/lss_mock_dev/python/esass/logs
source /home/erosita/sw/eSASSusers_200302/bin/esass-init.sh
nohup sh calibrate_events_0.sh > calibrate_events_batch0.log & 
nohup sh calibrate_events_1.sh > calibrate_events_batch1.log & 
nohup sh calibrate_events_2.sh > calibrate_events_batch2.log & 
nohup sh calibrate_events_3.sh > calibrate_events_batch3.log & 
nohup sh calibrate_events_4.sh > calibrate_events_batch4.log & 
nohup sh calibrate_events_5.sh > calibrate_events_batch5.log & 
nohup sh calibrate_events_6.sh > calibrate_events_batch6.log & 
nohup sh calibrate_events_7.sh > calibrate_events_batch7.log & 

/data26s/comparat/simulations/erosim/SFC-v1-erass1
/data26s/comparat/simulations/erosim/SFC-v1-erass8


python $GIT_AGN_MOCK/python/esass/calibrate_events.py 104 101.24999999999999 35.68533471265204


python $GIT_AGN_MOCK/python/esass/calibrate_events.py 172 199.28571428571428 48.141207794360284
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 185 213.75 66.44353569089877
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 186 191.25 66.44353569089877
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 187 195.0 72.38756092964962
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 188 224.99999999999997 72.38756092964962
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 191 225.0 84.14973293629666
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 203 303.75 24.62431835216408
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 216 326.24999999999994 35.68533471265204
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 306 354.375 9.594068226860458
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 307 0.0 14.477512185929925
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 309 16.875 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 310 5.625 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 311 11.25 24.62431835216408
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 314 343.12499999999994 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 317 5.625 29.999999999999993
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 319 0.0 35.68533471265204
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 320 90.0 -35.68533471265205
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 322 84.375 -30.000000000000014
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 323 90.0 -24.62431835216408
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 324 101.24999999999999 -24.62431835216408
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 325 106.875 -19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 326 95.625 -19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 327 101.24999999999999 -14.477512185929939
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 461 275.625 -9.594068226860458
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 497 275.625 9.594068226860458
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 500 281.24999999999994 14.477512185929925
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 501 286.875 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 502 275.625 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 503 281.24999999999994 24.62431835216408
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 504 258.75 14.477512185929925
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 505 264.37499999999994 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 506 253.12499999999994 19.47122063449069
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 507 258.75 24.62431835216408
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 509 275.625 29.999999999999993
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 510 264.37499999999994 29.999999999999993
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 511 270.0 35.68533471265204
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 512 45.0 -84.14973293629666
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 514 22.5 -78.28414760510762
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 515 45.0 -72.38756092964962
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 528 81.0 -60.4344388449523
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 532 83.57142857142858 -48.1412077943603
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 560 45.0 -35.68533471265205
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 604 157.5 -35.68533471265205


157
180
181
182
183
238
250
354
383
387
517
586
675
684
689
691
704
713
721
724
735


source /home/erosita/sw/sass-setup.sh eSASSusers_200302

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
ll /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000
ll /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528
ll /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608
ll /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699
ll /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709
ll /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 000 45.0 4.780191847199163
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 528 81.0 -60.4344388449523
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 608 98.99999999999999 -60.4344388449523
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 699 213.75 -14.477512185929939
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 709 348.75 -66.44353569089876
python $GIT_AGN_MOCK/python/esass/calibrate_events.py 767 315.0 -4.780191847199163


# The actual esass run 

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
cd $GIT_AGN_MOCK/python/esass
$GIT_AGN_MOCK/python/esass/eRO_pipe_images.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/evt_000.fits /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000 -b 16 --hp8
$GIT_AGN_MOCK/python/esass/eRO_pipe_det.py  /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000 -V 13

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
cd $GIT_AGN_MOCK/python/esass
$GIT_AGN_MOCK/python/esass/eRO_pipe_images.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699/evt_699.fits /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699 -b 16 --hp8
$GIT_AGN_MOCK/python/esass/eRO_pipe_det.py  /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699 -V 13

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
cd $GIT_AGN_MOCK/python/esass
$GIT_AGN_MOCK/python/esass/eRO_pipe_images.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767/evt_767.fits /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767 -b 16 --hp8
$GIT_AGN_MOCK/python/esass/eRO_pipe_det.py  /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767 -V 13

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
cd $GIT_AGN_MOCK/python/esass
$GIT_AGN_MOCK/python/esass/eRO_pipe_images.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709/evt_709.fits /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709 -b 16 --hp8
$GIT_AGN_MOCK/python/esass/eRO_pipe_det.py  /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709 -V 13

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
cd $GIT_AGN_MOCK/python/esass
$GIT_AGN_MOCK/python/esass/eRO_pipe_images.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608/evt_608.fits /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608 -b 16 --hp8
$GIT_AGN_MOCK/python/esass/eRO_pipe_det.py  /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608 -V 13

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
cd $GIT_AGN_MOCK/python/esass
$GIT_AGN_MOCK/python/esass/eRO_pipe_images.py /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528/evt_528.fits /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528 -b 16 --hp8
$GIT_AGN_MOCK/python/esass/eRO_pipe_det.py  /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528 /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528 -V 13

During the first command :

expmap/get_badpix: **ERROR2** WITHCALBADPIX set to false. Reading bad pixels from event file.

 IXY_SIZE_MERGED:       12000       12000
 expmap: DONE
 expmap: eSASSusers_190925 Okt 01 17:49:16 2019
 expmap/get_badpix: **ERROR2** WITHCALBADPIX set to false. Reading bad pixels from event file.
 expmap/ftiopn: **STOP** could not open the named file
 expmap: FAILED
Traceback (most recent call last):
  File "/home/comparat/software/linux/lss_mock_dev/python/esass/eRO_pipe_images.py", line 156, in <module>
    subprocess.run(cmd_expmap0, check=True)
  File "/home/comparat/miniconda3/envs/astroconda/lib/python3.7/subprocess.py", line 487, in run
    output=stdout, stderr=stderr)
subprocess.CalledProcessError: Command '['expmap', 'inputdatasets=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/evt_000.fits', 'emin=0.6', 'emax=2.3', 'templateimage=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/020_EvtImg.fits', 'withdetmaps=yes', 'withvignetting=no', 'withmergedmaps=yes', 'withfilebadpix=yes', 'withcalbadpix=yes', 'withweights=yes', 'mergedmaps=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/020_ExpMap.fits']' returned non-zero exit status 1.


# write the commands run eSASS
python run_esass_WriteCommands.py
# it creates 
# a set of 768 eRASS_det_pipe_V12.py commands written in 
esass_commands.sh
# run all the commands inside this script after having sourced eSASS
# it needs evtool eSASS command (i.e. to source eSASS)
source /home/erosita/sw/sass-setup.sh eSASSusers_190925

# Fields to start with for the SFC :
# 
grep 000 esass_commands.sh
grep 699 esass_commands.sh
grep 767 esass_commands.sh
grep 709 esass_commands.sh
grep 608 esass_commands.sh
grep 528 esass_commands.sh

source /home/erosita/sw/sass-setup.sh eSASSusers_190925
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/evt_000.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_v12_000_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699/evt_699.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_v12_699_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767/evt_767.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_v12_767_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709/evt_709.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_v12_709_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608/evt_608.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_v12_608_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528/evt_528.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_v12_528_eSASSusers_190925.log & 

nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12_final.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/evt_000.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_final_v12_000_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12_final.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699/evt_699.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_final_v12_699_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12_final.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767/evt_767.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_final_v12_767_eSASSusers_190925.log & 

nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12_final.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709/evt_709.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_final_v12_709_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12_final.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608/evt_608.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_final_v12_608_eSASSusers_190925.log & 
nohup python $GIT_AGN_MOCK/python/esass/eRASS_det_pipe_V12_final.py --infile=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528/evt_528.fits --outdir=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528/ --band=c > /home/comparat/software/linux/lss_mock_dev/python/esass/logs/esass_final_v12_528_eSASSusers_190925.log & 


# Rsync results to the SFC repository
# the catalogs
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/000 /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/528 /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/608 /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/699 /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/709 /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/767 /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/
# the photons
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/evt_???.fits /data42s/comparat/firefly/mocks/2020-03/SFC/events/

rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/0??_0216_EvtImg.fits    /data42s/comparat/firefly/mocks/2020-03/SFC/0216_EvtImg/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/0??_0216_UnvExpMap.fits /data42s/comparat/firefly/mocks/2020-03/SFC/0216_UnvExpMap/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/0??_0216_DetMsk.fits    /data42s/comparat/firefly/mocks/2020-03/SFC/0216_DetMsk/
rsync -avz /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/0??_0216_ExpMap.fits    /data42s/comparat/firefly/mocks/2020-03/SFC/0216_ExpMap/
# mkdir /data42s/comparat/firefly/mocks/2020-03/SFC/0216_EvtImg/
# mkdir /data42s/comparat/firefly/mocks/2020-03/SFC/0216_UnvExpMap/
# mkdir /data42s/comparat/firefly/mocks/2020-03/SFC/0216_DetMsk/
# mkdir /data42s/comparat/firefly/mocks/2020-03/SFC/0216_ExpMap/

# change permissions
chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/*
chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/events/*
chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/*
chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/eSASS_run/???/*

# rsync input catalogues
rsync -avz  $MD10/cat_AGN-MAG_all/*000000* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN-MAG_all/
rsync -avz  $MD10/cat_AGN-MAG_all/*000528* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN-MAG_all/
rsync -avz  $MD10/cat_AGN-MAG_all/*000608* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN-MAG_all/
rsync -avz  $MD10/cat_AGN-MAG_all/*000699* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN-MAG_all/
rsync -avz  $MD10/cat_AGN-MAG_all/*000709* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN-MAG_all/
rsync -avz  $MD10/cat_AGN-MAG_all/*000767* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN-MAG_all/

rsync -avz  $MD10/cat_AGN_SIMPUT/*000000* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN_SIMPUT/
rsync -avz  $MD10/cat_AGN_SIMPUT/*000528* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN_SIMPUT/
rsync -avz  $MD10/cat_AGN_SIMPUT/*000608* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN_SIMPUT/
rsync -avz  $MD10/cat_AGN_SIMPUT/*000699* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN_SIMPUT/
rsync -avz  $MD10/cat_AGN_SIMPUT/*000709* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN_SIMPUT/
rsync -avz  $MD10/cat_AGN_SIMPUT/*000767* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_AGN_SIMPUT/

rsync -avz  $MD10/cat_CLU_SIMPUT/*000000* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_CLU_SIMPUT/
rsync -avz  $MD10/cat_CLU_SIMPUT/*000528* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_CLU_SIMPUT/
rsync -avz  $MD10/cat_CLU_SIMPUT/*000608* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_CLU_SIMPUT/
rsync -avz  $MD10/cat_CLU_SIMPUT/*000699* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_CLU_SIMPUT/
rsync -avz  $MD10/cat_CLU_SIMPUT/*000709* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_CLU_SIMPUT/
rsync -avz  $MD10/cat_CLU_SIMPUT/*000767* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_CLU_SIMPUT/

rsync -avz  $MD10/cat_eRO_CLU/*000000* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_eRO_CLU/
rsync -avz  $MD10/cat_eRO_CLU/*000528* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_eRO_CLU/
rsync -avz  $MD10/cat_eRO_CLU/*000608* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_eRO_CLU/
rsync -avz  $MD10/cat_eRO_CLU/*000699* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_eRO_CLU/
rsync -avz  $MD10/cat_eRO_CLU/*000709* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_eRO_CLU/
rsync -avz  $MD10/cat_eRO_CLU/*000767* /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/cat_eRO_CLU/

chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/*
chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/*
chmod ugo+rX /data42s/comparat/firefly/mocks/2020-03/SFC/input_catalogues/*/*
