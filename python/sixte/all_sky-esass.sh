!/bin/bash 

on ds23

cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn
pyCONDA
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 

nohup sh cal_evt_000.sh > cal_evt_000.log &
nohup sh cal_evt_111.sh > cal_evt_111.log &
nohup sh cal_evt_222.sh > cal_evt_222.log &
nohup sh cal_evt_efeds.sh > cal_evt_efeds_c.log &

import numpy as n
for ii in n.arange(768):
	print("cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/"+str(ii).zfill(3))
	print("nohup ../eRASS_det_pipe_V12.py evt_"+str(ii).zfill(3)+".fits --bands=c > v12.log & ")

cd $GIT_AGN_MOCK/python/sixte/
ls /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/???/PipeV11.c_results/020_SrcCat.fits > $GIT_AGN_MOCK/python/sixte/srcCat.list

stilts tcat in=@srcCat.list ifmt=fits omode=out ofmt=fits out=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/merge_020_SrcCat.fits

stilts tcat in=@inputCat.list ifmt=fits omode=out ofmt=fits out=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/input_clusters.fits

stilts tcat in=@inputCatAGN.list ifmt=fits omode=out ofmt=fits out=/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/input_AGN.fits

ls /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/???/evt_???.fits > $GIT_AGN_MOCK/python/sixte/evt.list

ls ???/PipeV11.c_results/???_out_det-in.fits > out_det-in.list
stilts tcat in=@out_det-in.list ifmt=fits omode=out ofmt=fits out=out_det-in.fits

ls ???/PipeV11.c_results/???_out_in-det.fits > out_in-det.list
stilts tcat in=@out_in-det.list ifmt=fits omode=out ofmt=fits out=out_in-det.fits

# rm /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/???/PipeV11.c_results/???_out_det-in.fits
# rm /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/???/PipeV11.c_results/???_out_in-det.fits

cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/000
nohup ../eRASS_det_pipe_V12.py evt_000.fits --bands=c > v12_000.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/001     
# nohup ../eRASS_det_pipe_V12.py evt_001.fits --bands=c > v12_001.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/002     
# nohup ../eRASS_det_pipe_V12.py evt_002.fits --bands=c > v12_002.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/003     
# nohup ../eRASS_det_pipe_V12.py evt_003.fits --bands=c > v12_003.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/004     
# nohup ../eRASS_det_pipe_V12.py evt_004.fits --bands=c > v12_004.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/005     
# nohup ../eRASS_det_pipe_V12.py evt_005.fits --bands=c > v12_005.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/006     
# nohup ../eRASS_det_pipe_V12.py evt_006.fits --bands=c > v12_006.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/007     
# nohup ../eRASS_det_pipe_V12.py evt_007.fits --bands=c > v12_007.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/008     
# nohup ../eRASS_det_pipe_V12.py evt_008.fits --bands=c > v12_008.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/009     
# nohup ../eRASS_det_pipe_V12.py evt_009.fits --bands=c > v12_009.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/010     
# nohup ../eRASS_det_pipe_V12.py evt_010.fits --bands=c > v12_010.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/011     
# nohup ../eRASS_det_pipe_V12.py evt_011.fits --bands=c > v12_011.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/012     
# nohup ../eRASS_det_pipe_V12.py evt_012.fits --bands=c > v12_012.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/013     
# nohup ../eRASS_det_pipe_V12.py evt_013.fits --bands=c > v12_013.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/014     
# nohup ../eRASS_det_pipe_V12.py evt_014.fits --bands=c > v12_014.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/015     
# nohup ../eRASS_det_pipe_V12.py evt_015.fits --bands=c > v12_015.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/016     
# nohup ../eRASS_det_pipe_V12.py evt_016.fits --bands=c > v12_016.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/017     
# nohup ../eRASS_det_pipe_V12.py evt_017.fits --bands=c > v12_017.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/018     
# nohup ../eRASS_det_pipe_V12.py evt_018.fits --bands=c > v12_018.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/019     
# nohup ../eRASS_det_pipe_V12.py evt_019.fits --bands=c > v12_019.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/020     
# nohup ../eRASS_det_pipe_V12.py evt_020.fits --bands=c > v12_020.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/021     
# nohup ../eRASS_det_pipe_V12.py evt_021.fits --bands=c > v12_021.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/022     
# nohup ../eRASS_det_pipe_V12.py evt_022.fits --bands=c > v12_022.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/023     
# nohup ../eRASS_det_pipe_V12.py evt_023.fits --bands=c > v12_023.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/024     
# nohup ../eRASS_det_pipe_V12.py evt_024.fits --bands=c > v12_024.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/025     
# nohup ../eRASS_det_pipe_V12.py evt_025.fits --bands=c > v12_025.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/026     
# nohup ../eRASS_det_pipe_V12.py evt_026.fits --bands=c > v12_026.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/027     
# nohup ../eRASS_det_pipe_V12.py evt_027.fits --bands=c > v12_027.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/028     
# nohup ../eRASS_det_pipe_V12.py evt_028.fits --bands=c > v12_028.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/029     
# nohup ../eRASS_det_pipe_V12.py evt_029.fits --bands=c > v12_029.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/030     
# nohup ../eRASS_det_pipe_V12.py evt_030.fits --bands=c > v12_030.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/031     
# nohup ../eRASS_det_pipe_V12.py evt_031.fits --bands=c > v12_031.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/032     
# nohup ../eRASS_det_pipe_V12.py evt_032.fits --bands=c > v12_032.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/033     
# nohup ../eRASS_det_pipe_V12.py evt_033.fits --bands=c > v12_033.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/034     
# nohup ../eRASS_det_pipe_V12.py evt_034.fits --bands=c > v12_034.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/035     
# nohup ../eRASS_det_pipe_V12.py evt_035.fits --bands=c > v12_035.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/036     
# nohup ../eRASS_det_pipe_V12.py evt_036.fits --bands=c > v12_036.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/037     
# nohup ../eRASS_det_pipe_V12.py evt_037.fits --bands=c > v12_037.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/038     
# nohup ../eRASS_det_pipe_V12.py evt_038.fits --bands=c > v12_038.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/039     
# nohup ../eRASS_det_pipe_V12.py evt_039.fits --bands=c > v12_039.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/040     
# nohup ../eRASS_det_pipe_V12.py evt_040.fits --bands=c > v12_040.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/041     
# nohup ../eRASS_det_pipe_V12.py evt_041.fits --bands=c > v12_041.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/042     
# nohup ../eRASS_det_pipe_V12.py evt_042.fits --bands=c > v12_042.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/043     
# nohup ../eRASS_det_pipe_V12.py evt_043.fits --bands=c > v12_043.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/044     
# nohup ../eRASS_det_pipe_V12.py evt_044.fits --bands=c > v12_044.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/045     
# nohup ../eRASS_det_pipe_V12.py evt_045.fits --bands=c > v12_045.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/046     
# nohup ../eRASS_det_pipe_V12.py evt_046.fits --bands=c > v12_046.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/047     
# nohup ../eRASS_det_pipe_V12.py evt_047.fits --bands=c > v12_047.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/048     
# nohup ../eRASS_det_pipe_V12.py evt_048.fits --bands=c > v12_048.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/049     
# nohup ../eRASS_det_pipe_V12.py evt_049.fits --bands=c > v12_049.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/050     
# nohup ../eRASS_det_pipe_V12.py evt_050.fits --bands=c > v12_050.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/051     
# nohup ../eRASS_det_pipe_V12.py evt_051.fits --bands=c > v12_051.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/052     
# nohup ../eRASS_det_pipe_V12.py evt_052.fits --bands=c > v12_052.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/053     
# nohup ../eRASS_det_pipe_V12.py evt_053.fits --bands=c > v12_053.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/054     
# nohup ../eRASS_det_pipe_V12.py evt_054.fits --bands=c > v12_054.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/055     
# nohup ../eRASS_det_pipe_V12.py evt_055.fits --bands=c > v12_055.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/056     
# nohup ../eRASS_det_pipe_V12.py evt_056.fits --bands=c > v12_056.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/057     
# nohup ../eRASS_det_pipe_V12.py evt_057.fits --bands=c > v12_057.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/058     
# nohup ../eRASS_det_pipe_V12.py evt_058.fits --bands=c > v12_058.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/059     
# nohup ../eRASS_det_pipe_V12.py evt_059.fits --bands=c > v12_059.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/060     
# nohup ../eRASS_det_pipe_V12.py evt_060.fits --bands=c > v12_060.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/061     
# nohup ../eRASS_det_pipe_V12.py evt_061.fits --bands=c > v12_061.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/062     
# nohup ../eRASS_det_pipe_V12.py evt_062.fits --bands=c > v12_062.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/063     
# nohup ../eRASS_det_pipe_V12.py evt_063.fits --bands=c > v12_063.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/064     
# nohup ../eRASS_det_pipe_V12.py evt_064.fits --bands=c > v12_064.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/065     
# nohup ../eRASS_det_pipe_V12.py evt_065.fits --bands=c > v12_065.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/066     
# nohup ../eRASS_det_pipe_V12.py evt_066.fits --bands=c > v12_066.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/067     
# nohup ../eRASS_det_pipe_V12.py evt_067.fits --bands=c > v12_067.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/068     
# nohup ../eRASS_det_pipe_V12.py evt_068.fits --bands=c > v12_068.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/069     
# nohup ../eRASS_det_pipe_V12.py evt_069.fits --bands=c > v12_069.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/070     
# nohup ../eRASS_det_pipe_V12.py evt_070.fits --bands=c > v12_070.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/071     
# nohup ../eRASS_det_pipe_V12.py evt_071.fits --bands=c > v12_071.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/072     
# nohup ../eRASS_det_pipe_V12.py evt_072.fits --bands=c > v12_072.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/073     
# nohup ../eRASS_det_pipe_V12.py evt_073.fits --bands=c > v12_073.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/074     
# nohup ../eRASS_det_pipe_V12.py evt_074.fits --bands=c > v12_074.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/075     
# nohup ../eRASS_det_pipe_V12.py evt_075.fits --bands=c > v12_075.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/076     
# nohup ../eRASS_det_pipe_V12.py evt_076.fits --bands=c > v12_076.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/077     
# nohup ../eRASS_det_pipe_V12.py evt_077.fits --bands=c > v12_077.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/078     
# nohup ../eRASS_det_pipe_V12.py evt_078.fits --bands=c > v12_078.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/079     
# nohup ../eRASS_det_pipe_V12.py evt_079.fits --bands=c > v12_079.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/080     
# nohup ../eRASS_det_pipe_V12.py evt_080.fits --bands=c > v12_080.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/081     
# nohup ../eRASS_det_pipe_V12.py evt_081.fits --bands=c > v12_081.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/082     
# nohup ../eRASS_det_pipe_V12.py evt_082.fits --bands=c > v12_082.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/083     
# nohup ../eRASS_det_pipe_V12.py evt_083.fits --bands=c > v12_083.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/084     
# nohup ../eRASS_det_pipe_V12.py evt_084.fits --bands=c > v12_084.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/085     
# nohup ../eRASS_det_pipe_V12.py evt_085.fits --bands=c > v12_085.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/086     
# nohup ../eRASS_det_pipe_V12.py evt_086.fits --bands=c > v12_086.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/087     
# nohup ../eRASS_det_pipe_V12.py evt_087.fits --bands=c > v12_087.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/088     
# nohup ../eRASS_det_pipe_V12.py evt_088.fits --bands=c > v12_088.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/089     
# nohup ../eRASS_det_pipe_V12.py evt_089.fits --bands=c > v12_089.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/090     
# nohup ../eRASS_det_pipe_V12.py evt_090.fits --bands=c > v12_090.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/091     
# nohup ../eRASS_det_pipe_V12.py evt_091.fits --bands=c > v12_091.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/092     
# nohup ../eRASS_det_pipe_V12.py evt_092.fits --bands=c > v12_092.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/093     
# nohup ../eRASS_det_pipe_V12.py evt_093.fits --bands=c > v12_093.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/094     
# nohup ../eRASS_det_pipe_V12.py evt_094.fits --bands=c > v12_094.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/095     
# nohup ../eRASS_det_pipe_V12.py evt_095.fits --bands=c > v12_095.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/096     
# nohup ../eRASS_det_pipe_V12.py evt_096.fits --bands=c > v12_096.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/097     
# nohup ../eRASS_det_pipe_V12.py evt_097.fits --bands=c > v12_097.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/098     
# nohup ../eRASS_det_pipe_V12.py evt_098.fits --bands=c > v12_098.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/099     
# nohup ../eRASS_det_pipe_V12.py evt_099.fits --bands=c > v12_099.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/100     
# nohup ../eRASS_det_pipe_V12.py evt_100.fits --bands=c > v12_100.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/101     
# nohup ../eRASS_det_pipe_V12.py evt_101.fits --bands=c > v12_101.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/102     
# nohup ../eRASS_det_pipe_V12.py evt_102.fits --bands=c > v12_102.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/103     
# nohup ../eRASS_det_pipe_V12.py evt_103.fits --bands=c > v12_103.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/104     
# nohup ../eRASS_det_pipe_V12.py evt_104.fits --bands=c > v12_104.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/105     
# nohup ../eRASS_det_pipe_V12.py evt_105.fits --bands=c > v12_105.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/106     
# nohup ../eRASS_det_pipe_V12.py evt_106.fits --bands=c > v12_106.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/107     
# nohup ../eRASS_det_pipe_V12.py evt_107.fits --bands=c > v12_107.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/108     
# nohup ../eRASS_det_pipe_V12.py evt_108.fits --bands=c > v12_108.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/109     
# nohup ../eRASS_det_pipe_V12.py evt_109.fits --bands=c > v12_109.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/110     
# nohup ../eRASS_det_pipe_V12.py evt_110.fits --bands=c > v12_110.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/111     
# nohup ../eRASS_det_pipe_V12.py evt_111.fits --bands=c > v12_111.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/112     
# nohup ../eRASS_det_pipe_V12.py evt_112.fits --bands=c > v12_112.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/113     
# nohup ../eRASS_det_pipe_V12.py evt_113.fits --bands=c > v12_113.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/114     
# nohup ../eRASS_det_pipe_V12.py evt_114.fits --bands=c > v12_114.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/115     
# nohup ../eRASS_det_pipe_V12.py evt_115.fits --bands=c > v12_115.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/116     
# nohup ../eRASS_det_pipe_V12.py evt_116.fits --bands=c > v12_116.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/117     
# nohup ../eRASS_det_pipe_V12.py evt_117.fits --bands=c > v12_117.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/118     
# nohup ../eRASS_det_pipe_V12.py evt_118.fits --bands=c > v12_118.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/119     
# nohup ../eRASS_det_pipe_V12.py evt_119.fits --bands=c > v12_119.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/120     
# nohup ../eRASS_det_pipe_V12.py evt_120.fits --bands=c > v12_120.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/121     
# nohup ../eRASS_det_pipe_V12.py evt_121.fits --bands=c > v12_121.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/122     
# nohup ../eRASS_det_pipe_V12.py evt_122.fits --bands=c > v12_122.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/123     
# nohup ../eRASS_det_pipe_V12.py evt_123.fits --bands=c > v12_123.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/124     
# nohup ../eRASS_det_pipe_V12.py evt_124.fits --bands=c > v12_124.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/125     
# nohup ../eRASS_det_pipe_V12.py evt_125.fits --bands=c > v12_125.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/126     
# nohup ../eRASS_det_pipe_V12.py evt_126.fits --bands=c > v12_126.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/127     
# nohup ../eRASS_det_pipe_V12.py evt_127.fits --bands=c > v12_127.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/128     
# nohup ../eRASS_det_pipe_V12.py evt_128.fits --bands=c > v12_128.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/129     
# nohup ../eRASS_det_pipe_V12.py evt_129.fits --bands=c > v12_129.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/130     
# nohup ../eRASS_det_pipe_V12.py evt_130.fits --bands=c > v12_130.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/131     
# nohup ../eRASS_det_pipe_V12.py evt_131.fits --bands=c > v12_131.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/132     
# nohup ../eRASS_det_pipe_V12.py evt_132.fits --bands=c > v12_132.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/133     
# nohup ../eRASS_det_pipe_V12.py evt_133.fits --bands=c > v12_133.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/134     
# nohup ../eRASS_det_pipe_V12.py evt_134.fits --bands=c > v12_134.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/135     
# nohup ../eRASS_det_pipe_V12.py evt_135.fits --bands=c > v12_135.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/136     
# nohup ../eRASS_det_pipe_V12.py evt_136.fits --bands=c > v12_136.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/137     
# nohup ../eRASS_det_pipe_V12.py evt_137.fits --bands=c > v12_137.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/138     
# nohup ../eRASS_det_pipe_V12.py evt_138.fits --bands=c > v12_138.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/139     
# nohup ../eRASS_det_pipe_V12.py evt_139.fits --bands=c > v12_139.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/140     
# nohup ../eRASS_det_pipe_V12.py evt_140.fits --bands=c > v12_140.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/141     
# nohup ../eRASS_det_pipe_V12.py evt_141.fits --bands=c > v12_141.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/142     
# nohup ../eRASS_det_pipe_V12.py evt_142.fits --bands=c > v12_142.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/143     
# nohup ../eRASS_det_pipe_V12.py evt_143.fits --bands=c > v12_143.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/144     
# nohup ../eRASS_det_pipe_V12.py evt_144.fits --bands=c > v12_144.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/145     
# nohup ../eRASS_det_pipe_V12.py evt_145.fits --bands=c > v12_145.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/146     
# nohup ../eRASS_det_pipe_V12.py evt_146.fits --bands=c > v12_146.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/147     
# nohup ../eRASS_det_pipe_V12.py evt_147.fits --bands=c > v12_147.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/148     
# nohup ../eRASS_det_pipe_V12.py evt_148.fits --bands=c > v12_148.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/149     
# nohup ../eRASS_det_pipe_V12.py evt_149.fits --bands=c > v12_149.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/150     
# nohup ../eRASS_det_pipe_V12.py evt_150.fits --bands=c > v12_150.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/151     
nohup ../eRASS_det_pipe_V12.py evt_151.fits --bands=c > v12_151.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/152     
nohup ../eRASS_det_pipe_V12.py evt_152.fits --bands=c > v12_152.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/153     
nohup ../eRASS_det_pipe_V12.py evt_153.fits --bands=c > v12_153.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/154     
nohup ../eRASS_det_pipe_V12.py evt_154.fits --bands=c > v12_154.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/155     
nohup ../eRASS_det_pipe_V12.py evt_155.fits --bands=c > v12_155.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/156     
nohup ../eRASS_det_pipe_V12.py evt_156.fits --bands=c > v12_156.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/157     
nohup ../eRASS_det_pipe_V12.py evt_157.fits --bands=c > v12_157.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/158     
nohup ../eRASS_det_pipe_V12.py evt_158.fits --bands=c > v12_158.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/159     
nohup ../eRASS_det_pipe_V12.py evt_159.fits --bands=c > v12_159.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/160     
# nohup ../eRASS_det_pipe_V12.py evt_160.fits --bands=c > v12_160.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/161     
# nohup ../eRASS_det_pipe_V12.py evt_161.fits --bands=c > v12_161.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/162     
# nohup ../eRASS_det_pipe_V12.py evt_162.fits --bands=c > v12_162.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/163     
# nohup ../eRASS_det_pipe_V12.py evt_163.fits --bands=c > v12_163.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/164     
# nohup ../eRASS_det_pipe_V12.py evt_164.fits --bands=c > v12_164.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/165     
# nohup ../eRASS_det_pipe_V12.py evt_165.fits --bands=c > v12_165.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/166     
# nohup ../eRASS_det_pipe_V12.py evt_166.fits --bands=c > v12_166.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/167     
# nohup ../eRASS_det_pipe_V12.py evt_167.fits --bands=c > v12_167.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/168     
# nohup ../eRASS_det_pipe_V12.py evt_168.fits --bands=c > v12_168.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/169     
# nohup ../eRASS_det_pipe_V12.py evt_169.fits --bands=c > v12_169.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/170     
# nohup ../eRASS_det_pipe_V12.py evt_170.fits --bands=c > v12_170.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/171     
# nohup ../eRASS_det_pipe_V12.py evt_171.fits --bands=c > v12_171.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/172     
# nohup ../eRASS_det_pipe_V12.py evt_172.fits --bands=c > v12_172.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/173     
# nohup ../eRASS_det_pipe_V12.py evt_173.fits --bands=c > v12_173.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/174     
# nohup ../eRASS_det_pipe_V12.py evt_174.fits --bands=c > v12_174.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/175     
# nohup ../eRASS_det_pipe_V12.py evt_175.fits --bands=c > v12_175.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/176     
# nohup ../eRASS_det_pipe_V12.py evt_176.fits --bands=c > v12_176.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/177     
# nohup ../eRASS_det_pipe_V12.py evt_177.fits --bands=c > v12_177.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/178     
# nohup ../eRASS_det_pipe_V12.py evt_178.fits --bands=c > v12_178.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/179     
# nohup ../eRASS_det_pipe_V12.py evt_179.fits --bands=c > v12_179.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/180     
nohup ../eRASS_det_pipe_V12.py evt_180.fits --bands=c > v12_180.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/181     
nohup ../eRASS_det_pipe_V12.py evt_181.fits --bands=c > v12_181.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/182     
nohup ../eRASS_det_pipe_V12.py evt_182.fits --bands=c > v12_182.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/183     
nohup ../eRASS_det_pipe_V12.py evt_183.fits --bands=c > v12_183.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/184     
nohup ../eRASS_det_pipe_V12.py evt_184.fits --bands=c > v12_184.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/185     
# nohup ../eRASS_det_pipe_V12.py evt_185.fits --bands=c > v12_185.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/186     
# nohup ../eRASS_det_pipe_V12.py evt_186.fits --bands=c > v12_186.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/187     
# nohup ../eRASS_det_pipe_V12.py evt_187.fits --bands=c > v12_187.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/188     
# nohup ../eRASS_det_pipe_V12.py evt_188.fits --bands=c > v12_188.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/189     
# nohup ../eRASS_det_pipe_V12.py evt_189.fits --bands=c > v12_189.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/190     
# nohup ../eRASS_det_pipe_V12.py evt_190.fits --bands=c > v12_190.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/191     
# nohup ../eRASS_det_pipe_V12.py evt_191.fits --bands=c > v12_191.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/192     
# nohup ../eRASS_det_pipe_V12.py evt_192.fits --bands=c > v12_192.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/193     
# nohup ../eRASS_det_pipe_V12.py evt_193.fits --bands=c > v12_193.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/194     
# nohup ../eRASS_det_pipe_V12.py evt_194.fits --bands=c > v12_194.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/195     
# nohup ../eRASS_det_pipe_V12.py evt_195.fits --bands=c > v12_195.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/196     
# nohup ../eRASS_det_pipe_V12.py evt_196.fits --bands=c > v12_196.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/197     
# nohup ../eRASS_det_pipe_V12.py evt_197.fits --bands=c > v12_197.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/198     
# nohup ../eRASS_det_pipe_V12.py evt_198.fits --bands=c > v12_198.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/199     
# nohup ../eRASS_det_pipe_V12.py evt_199.fits --bands=c > v12_199.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/200     
# nohup ../eRASS_det_pipe_V12.py evt_200.fits --bands=c > v12_200.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/201     
# nohup ../eRASS_det_pipe_V12.py evt_201.fits --bands=c > v12_201.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/202     
# nohup ../eRASS_det_pipe_V12.py evt_202.fits --bands=c > v12_202.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/203     
# nohup ../eRASS_det_pipe_V12.py evt_203.fits --bands=c > v12_203.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/204     
# nohup ../eRASS_det_pipe_V12.py evt_204.fits --bands=c > v12_204.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/205     
# nohup ../eRASS_det_pipe_V12.py evt_205.fits --bands=c > v12_205.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/206     
# nohup ../eRASS_det_pipe_V12.py evt_206.fits --bands=c > v12_206.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/207     
# nohup ../eRASS_det_pipe_V12.py evt_207.fits --bands=c > v12_207.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/208     
# nohup ../eRASS_det_pipe_V12.py evt_208.fits --bands=c > v12_208.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/209     
# nohup ../eRASS_det_pipe_V12.py evt_209.fits --bands=c > v12_209.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/210     
# nohup ../eRASS_det_pipe_V12.py evt_210.fits --bands=c > v12_210.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/211     
# nohup ../eRASS_det_pipe_V12.py evt_211.fits --bands=c > v12_211.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/212     
# nohup ../eRASS_det_pipe_V12.py evt_212.fits --bands=c > v12_212.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/213     
# nohup ../eRASS_det_pipe_V12.py evt_213.fits --bands=c > v12_213.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/214     
# nohup ../eRASS_det_pipe_V12.py evt_214.fits --bands=c > v12_214.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/215     
# nohup ../eRASS_det_pipe_V12.py evt_215.fits --bands=c > v12_215.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/216     
# nohup ../eRASS_det_pipe_V12.py evt_216.fits --bands=c > v12_216.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/217     
# nohup ../eRASS_det_pipe_V12.py evt_217.fits --bands=c > v12_217.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/218     
# nohup ../eRASS_det_pipe_V12.py evt_218.fits --bands=c > v12_218.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/219     
# nohup ../eRASS_det_pipe_V12.py evt_219.fits --bands=c > v12_219.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/220     
# nohup ../eRASS_det_pipe_V12.py evt_220.fits --bands=c > v12_220.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/221     
# nohup ../eRASS_det_pipe_V12.py evt_221.fits --bands=c > v12_221.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/222     
# nohup ../eRASS_det_pipe_V12.py evt_222.fits --bands=c > v12_222.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/223     
# nohup ../eRASS_det_pipe_V12.py evt_223.fits --bands=c > v12_223.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/224     
# nohup ../eRASS_det_pipe_V12.py evt_224.fits --bands=c > v12_224.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/225     
# nohup ../eRASS_det_pipe_V12.py evt_225.fits --bands=c > v12_225.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/226     
# nohup ../eRASS_det_pipe_V12.py evt_226.fits --bands=c > v12_226.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/227     
# nohup ../eRASS_det_pipe_V12.py evt_227.fits --bands=c > v12_227.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/228     
# nohup ../eRASS_det_pipe_V12.py evt_228.fits --bands=c > v12_228.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/229     
# nohup ../eRASS_det_pipe_V12.py evt_229.fits --bands=c > v12_229.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/230     
# nohup ../eRASS_det_pipe_V12.py evt_230.fits --bands=c > v12_230.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/231     
# nohup ../eRASS_det_pipe_V12.py evt_231.fits --bands=c > v12_231.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/232     
# nohup ../eRASS_det_pipe_V12.py evt_232.fits --bands=c > v12_232.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/233     
# nohup ../eRASS_det_pipe_V12.py evt_233.fits --bands=c > v12_233.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/234     
# nohup ../eRASS_det_pipe_V12.py evt_234.fits --bands=c > v12_234.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/235     
# nohup ../eRASS_det_pipe_V12.py evt_235.fits --bands=c > v12_235.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/236     
# nohup ../eRASS_det_pipe_V12.py evt_236.fits --bands=c > v12_236.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/237     
# nohup ../eRASS_det_pipe_V12.py evt_237.fits --bands=c > v12_237.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/238     
nohup ../eRASS_det_pipe_V12.py evt_238.fits --bands=c > v12_238.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/239     
nohup ../eRASS_det_pipe_V12.py evt_239.fits --bands=c > v12_239.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/240     
# nohup ../eRASS_det_pipe_V12.py evt_240.fits --bands=c > v12_240.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/241     
# nohup ../eRASS_det_pipe_V12.py evt_241.fits --bands=c > v12_241.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/242     
# nohup ../eRASS_det_pipe_V12.py evt_242.fits --bands=c > v12_242.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/243     
# nohup ../eRASS_det_pipe_V12.py evt_243.fits --bands=c > v12_243.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/244     
# nohup ../eRASS_det_pipe_V12.py evt_244.fits --bands=c > v12_244.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/245     
# nohup ../eRASS_det_pipe_V12.py evt_245.fits --bands=c > v12_245.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/246     
# nohup ../eRASS_det_pipe_V12.py evt_246.fits --bands=c > v12_246.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/247     
# nohup ../eRASS_det_pipe_V12.py evt_247.fits --bands=c > v12_247.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/248     
# nohup ../eRASS_det_pipe_V12.py evt_248.fits --bands=c > v12_248.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/249     
# nohup ../eRASS_det_pipe_V12.py evt_249.fits --bands=c > v12_249.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/250     
# nohup ../eRASS_det_pipe_V12.py evt_250.fits --bands=c > v12_250.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/251     
# nohup ../eRASS_det_pipe_V12.py evt_251.fits --bands=c > v12_251.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/252     
# nohup ../eRASS_det_pipe_V12.py evt_252.fits --bands=c > v12_252.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/253     
# nohup ../eRASS_det_pipe_V12.py evt_253.fits --bands=c > v12_253.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/254     
# nohup ../eRASS_det_pipe_V12.py evt_254.fits --bands=c > v12_254.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/255     
# nohup ../eRASS_det_pipe_V12.py evt_255.fits --bands=c > v12_255.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/256     
# nohup ../eRASS_det_pipe_V12.py evt_256.fits --bands=c > v12_256.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/257     
# nohup ../eRASS_det_pipe_V12.py evt_257.fits --bands=c > v12_257.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/258     
# nohup ../eRASS_det_pipe_V12.py evt_258.fits --bands=c > v12_258.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/259     
# nohup ../eRASS_det_pipe_V12.py evt_259.fits --bands=c > v12_259.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/260     
# nohup ../eRASS_det_pipe_V12.py evt_260.fits --bands=c > v12_260.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/261     
# nohup ../eRASS_det_pipe_V12.py evt_261.fits --bands=c > v12_261.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/262     
# nohup ../eRASS_det_pipe_V12.py evt_262.fits --bands=c > v12_262.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/263     
# nohup ../eRASS_det_pipe_V12.py evt_263.fits --bands=c > v12_263.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/264     
# nohup ../eRASS_det_pipe_V12.py evt_264.fits --bands=c > v12_264.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/265     
# nohup ../eRASS_det_pipe_V12.py evt_265.fits --bands=c > v12_265.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/266     
# nohup ../eRASS_det_pipe_V12.py evt_266.fits --bands=c > v12_266.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/267     
# nohup ../eRASS_det_pipe_V12.py evt_267.fits --bands=c > v12_267.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/268     
# nohup ../eRASS_det_pipe_V12.py evt_268.fits --bands=c > v12_268.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/269     
# nohup ../eRASS_det_pipe_V12.py evt_269.fits --bands=c > v12_269.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/270     
# nohup ../eRASS_det_pipe_V12.py evt_270.fits --bands=c > v12_270.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/271     
# nohup ../eRASS_det_pipe_V12.py evt_271.fits --bands=c > v12_271.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/272     
# nohup ../eRASS_det_pipe_V12.py evt_272.fits --bands=c > v12_272.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/273     
# nohup ../eRASS_det_pipe_V12.py evt_273.fits --bands=c > v12_273.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/274     
# nohup ../eRASS_det_pipe_V12.py evt_274.fits --bands=c > v12_274.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/275     
# nohup ../eRASS_det_pipe_V12.py evt_275.fits --bands=c > v12_275.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/276     
# nohup ../eRASS_det_pipe_V12.py evt_276.fits --bands=c > v12_276.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/277     
# nohup ../eRASS_det_pipe_V12.py evt_277.fits --bands=c > v12_277.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/278     
# nohup ../eRASS_det_pipe_V12.py evt_278.fits --bands=c > v12_278.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/279     
# nohup ../eRASS_det_pipe_V12.py evt_279.fits --bands=c > v12_279.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/280     
# nohup ../eRASS_det_pipe_V12.py evt_280.fits --bands=c > v12_280.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/281     
# nohup ../eRASS_det_pipe_V12.py evt_281.fits --bands=c > v12_281.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/282     
# nohup ../eRASS_det_pipe_V12.py evt_282.fits --bands=c > v12_282.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/283     
# nohup ../eRASS_det_pipe_V12.py evt_283.fits --bands=c > v12_283.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/284     
# nohup ../eRASS_det_pipe_V12.py evt_284.fits --bands=c > v12_284.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/285     
# nohup ../eRASS_det_pipe_V12.py evt_285.fits --bands=c > v12_285.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/286     
# nohup ../eRASS_det_pipe_V12.py evt_286.fits --bands=c > v12_286.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/287     
# nohup ../eRASS_det_pipe_V12.py evt_287.fits --bands=c > v12_287.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/288     
# nohup ../eRASS_det_pipe_V12.py evt_288.fits --bands=c > v12_288.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/289     
# nohup ../eRASS_det_pipe_V12.py evt_289.fits --bands=c > v12_289.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/290     
# nohup ../eRASS_det_pipe_V12.py evt_290.fits --bands=c > v12_290.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/291     
# nohup ../eRASS_det_pipe_V12.py evt_291.fits --bands=c > v12_291.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/292     
# nohup ../eRASS_det_pipe_V12.py evt_292.fits --bands=c > v12_292.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/293     
# nohup ../eRASS_det_pipe_V12.py evt_293.fits --bands=c > v12_293.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/294     
# nohup ../eRASS_det_pipe_V12.py evt_294.fits --bands=c > v12_294.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/295     
# nohup ../eRASS_det_pipe_V12.py evt_295.fits --bands=c > v12_295.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/296     
# nohup ../eRASS_det_pipe_V12.py evt_296.fits --bands=c > v12_296.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/297     
# nohup ../eRASS_det_pipe_V12.py evt_297.fits --bands=c > v12_297.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/298     
# nohup ../eRASS_det_pipe_V12.py evt_298.fits --bands=c > v12_298.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/299     
# nohup ../eRASS_det_pipe_V12.py evt_299.fits --bands=c > v12_299.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/300     
# nohup ../eRASS_det_pipe_V12.py evt_300.fits --bands=c > v12_300.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/301     
# nohup ../eRASS_det_pipe_V12.py evt_301.fits --bands=c > v12_301.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/302     
# nohup ../eRASS_det_pipe_V12.py evt_302.fits --bands=c > v12_302.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/303     
# nohup ../eRASS_det_pipe_V12.py evt_303.fits --bands=c > v12_303.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/304     
# nohup ../eRASS_det_pipe_V12.py evt_304.fits --bands=c > v12_304.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/305     
# nohup ../eRASS_det_pipe_V12.py evt_305.fits --bands=c > v12_305.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/306     
# nohup ../eRASS_det_pipe_V12.py evt_306.fits --bands=c > v12_306.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/307     
# nohup ../eRASS_det_pipe_V12.py evt_307.fits --bands=c > v12_307.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/308     
# nohup ../eRASS_det_pipe_V12.py evt_308.fits --bands=c > v12_308.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/309     
# nohup ../eRASS_det_pipe_V12.py evt_309.fits --bands=c > v12_309.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/310     
# nohup ../eRASS_det_pipe_V12.py evt_310.fits --bands=c > v12_310.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/311     
# nohup ../eRASS_det_pipe_V12.py evt_311.fits --bands=c > v12_311.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/312     
# nohup ../eRASS_det_pipe_V12.py evt_312.fits --bands=c > v12_312.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/313     
# nohup ../eRASS_det_pipe_V12.py evt_313.fits --bands=c > v12_313.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/314     
# nohup ../eRASS_det_pipe_V12.py evt_314.fits --bands=c > v12_314.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/315     
# nohup ../eRASS_det_pipe_V12.py evt_315.fits --bands=c > v12_315.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/316     
# nohup ../eRASS_det_pipe_V12.py evt_316.fits --bands=c > v12_316.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/317     
# nohup ../eRASS_det_pipe_V12.py evt_317.fits --bands=c > v12_317.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/318     
# nohup ../eRASS_det_pipe_V12.py evt_318.fits --bands=c > v12_318.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/319     
# nohup ../eRASS_det_pipe_V12.py evt_319.fits --bands=c > v12_319.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/320     
# nohup ../eRASS_det_pipe_V12.py evt_320.fits --bands=c > v12_320.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/321     
# nohup ../eRASS_det_pipe_V12.py evt_321.fits --bands=c > v12_321.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/322     
# nohup ../eRASS_det_pipe_V12.py evt_322.fits --bands=c > v12_322.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/323     
# nohup ../eRASS_det_pipe_V12.py evt_323.fits --bands=c > v12_323.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/324     
# nohup ../eRASS_det_pipe_V12.py evt_324.fits --bands=c > v12_324.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/325     
# nohup ../eRASS_det_pipe_V12.py evt_325.fits --bands=c > v12_325.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/326     
# nohup ../eRASS_det_pipe_V12.py evt_326.fits --bands=c > v12_326.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/327     
# nohup ../eRASS_det_pipe_V12.py evt_327.fits --bands=c > v12_327.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/328     
# nohup ../eRASS_det_pipe_V12.py evt_328.fits --bands=c > v12_328.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/329     
# nohup ../eRASS_det_pipe_V12.py evt_329.fits --bands=c > v12_329.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/330     
# nohup ../eRASS_det_pipe_V12.py evt_330.fits --bands=c > v12_330.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/331     
# nohup ../eRASS_det_pipe_V12.py evt_331.fits --bands=c > v12_331.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/332     
# nohup ../eRASS_det_pipe_V12.py evt_332.fits --bands=c > v12_332.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/333     
# nohup ../eRASS_det_pipe_V12.py evt_333.fits --bands=c > v12_333.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/334     
# nohup ../eRASS_det_pipe_V12.py evt_334.fits --bands=c > v12_334.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/335     
# nohup ../eRASS_det_pipe_V12.py evt_335.fits --bands=c > v12_335.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/336     
# nohup ../eRASS_det_pipe_V12.py evt_336.fits --bands=c > v12_336.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/337     
# nohup ../eRASS_det_pipe_V12.py evt_337.fits --bands=c > v12_337.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/338     
# nohup ../eRASS_det_pipe_V12.py evt_338.fits --bands=c > v12_338.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/339     
# nohup ../eRASS_det_pipe_V12.py evt_339.fits --bands=c > v12_339.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/340     
# nohup ../eRASS_det_pipe_V12.py evt_340.fits --bands=c > v12_340.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/341     
# nohup ../eRASS_det_pipe_V12.py evt_341.fits --bands=c > v12_341.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/342     
# nohup ../eRASS_det_pipe_V12.py evt_342.fits --bands=c > v12_342.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/343     
# nohup ../eRASS_det_pipe_V12.py evt_343.fits --bands=c > v12_343.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/344     
# nohup ../eRASS_det_pipe_V12.py evt_344.fits --bands=c > v12_344.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/345     
# nohup ../eRASS_det_pipe_V12.py evt_345.fits --bands=c > v12_345.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/346     
# nohup ../eRASS_det_pipe_V12.py evt_346.fits --bands=c > v12_346.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/347     
# nohup ../eRASS_det_pipe_V12.py evt_347.fits --bands=c > v12_347.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/348     
# nohup ../eRASS_det_pipe_V12.py evt_348.fits --bands=c > v12_348.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/349     
# nohup ../eRASS_det_pipe_V12.py evt_349.fits --bands=c > v12_349.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/350     
# nohup ../eRASS_det_pipe_V12.py evt_350.fits --bands=c > v12_350.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/351     
# nohup ../eRASS_det_pipe_V12.py evt_351.fits --bands=c > v12_351.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/352     
# nohup ../eRASS_det_pipe_V12.py evt_352.fits --bands=c > v12_352.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/353     
# nohup ../eRASS_det_pipe_V12.py evt_353.fits --bands=c > v12_353.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/354     
# nohup ../eRASS_det_pipe_V12.py evt_354.fits --bands=c > v12_354.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/355     
# nohup ../eRASS_det_pipe_V12.py evt_355.fits --bands=c > v12_355.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/356     
# nohup ../eRASS_det_pipe_V12.py evt_356.fits --bands=c > v12_356.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/357     
# nohup ../eRASS_det_pipe_V12.py evt_357.fits --bands=c > v12_357.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/358     
# nohup ../eRASS_det_pipe_V12.py evt_358.fits --bands=c > v12_358.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/359     
# nohup ../eRASS_det_pipe_V12.py evt_359.fits --bands=c > v12_359.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/360     
# nohup ../eRASS_det_pipe_V12.py evt_360.fits --bands=c > v12_360.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/361     
# nohup ../eRASS_det_pipe_V12.py evt_361.fits --bands=c > v12_361.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/362     
# nohup ../eRASS_det_pipe_V12.py evt_362.fits --bands=c > v12_362.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/363     
# nohup ../eRASS_det_pipe_V12.py evt_363.fits --bands=c > v12_363.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/364     
# nohup ../eRASS_det_pipe_V12.py evt_364.fits --bands=c > v12_364.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/365     
# nohup ../eRASS_det_pipe_V12.py evt_365.fits --bands=c > v12_365.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/366     
# nohup ../eRASS_det_pipe_V12.py evt_366.fits --bands=c > v12_366.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/367     
# nohup ../eRASS_det_pipe_V12.py evt_367.fits --bands=c > v12_367.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/368     
# nohup ../eRASS_det_pipe_V12.py evt_368.fits --bands=c > v12_368.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/369     
# nohup ../eRASS_det_pipe_V12.py evt_369.fits --bands=c > v12_369.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/370     
# nohup ../eRASS_det_pipe_V12.py evt_370.fits --bands=c > v12_370.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/371     
# nohup ../eRASS_det_pipe_V12.py evt_371.fits --bands=c > v12_371.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/372     
# nohup ../eRASS_det_pipe_V12.py evt_372.fits --bands=c > v12_372.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/373     
# nohup ../eRASS_det_pipe_V12.py evt_373.fits --bands=c > v12_373.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/374     
# nohup ../eRASS_det_pipe_V12.py evt_374.fits --bands=c > v12_374.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/375     
# nohup ../eRASS_det_pipe_V12.py evt_375.fits --bands=c > v12_375.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/376     
# nohup ../eRASS_det_pipe_V12.py evt_376.fits --bands=c > v12_376.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/377     
# nohup ../eRASS_det_pipe_V12.py evt_377.fits --bands=c > v12_377.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/378     
# nohup ../eRASS_det_pipe_V12.py evt_378.fits --bands=c > v12_378.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/379     
# nohup ../eRASS_det_pipe_V12.py evt_379.fits --bands=c > v12_379.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/380     
# nohup ../eRASS_det_pipe_V12.py evt_380.fits --bands=c > v12_380.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/381     
# nohup ../eRASS_det_pipe_V12.py evt_381.fits --bands=c > v12_381.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/382     
# nohup ../eRASS_det_pipe_V12.py evt_382.fits --bands=c > v12_382.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/383     
# nohup ../eRASS_det_pipe_V12.py evt_383.fits --bands=c > v12_383.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/384     
# nohup ../eRASS_det_pipe_V12.py evt_384.fits --bands=c > v12_384.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/385     
# nohup ../eRASS_det_pipe_V12.py evt_385.fits --bands=c > v12_385.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/386     
# nohup ../eRASS_det_pipe_V12.py evt_386.fits --bands=c > v12_386.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/387     
# nohup ../eRASS_det_pipe_V12.py evt_387.fits --bands=c > v12_387.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/388     
# nohup ../eRASS_det_pipe_V12.py evt_388.fits --bands=c > v12_388.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/389     
# nohup ../eRASS_det_pipe_V12.py evt_389.fits --bands=c > v12_389.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/390     
# nohup ../eRASS_det_pipe_V12.py evt_390.fits --bands=c > v12_390.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/391     
# nohup ../eRASS_det_pipe_V12.py evt_391.fits --bands=c > v12_391.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/392     
# nohup ../eRASS_det_pipe_V12.py evt_392.fits --bands=c > v12_392.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/393     
# nohup ../eRASS_det_pipe_V12.py evt_393.fits --bands=c > v12_393.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/394     
# nohup ../eRASS_det_pipe_V12.py evt_394.fits --bands=c > v12_394.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/395     
# nohup ../eRASS_det_pipe_V12.py evt_395.fits --bands=c > v12_395.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/396     
# nohup ../eRASS_det_pipe_V12.py evt_396.fits --bands=c > v12_396.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/397     
# nohup ../eRASS_det_pipe_V12.py evt_397.fits --bands=c > v12_397.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/398     
# nohup ../eRASS_det_pipe_V12.py evt_398.fits --bands=c > v12_398.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/399     
# nohup ../eRASS_det_pipe_V12.py evt_399.fits --bands=c > v12_399.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/400     
# nohup ../eRASS_det_pipe_V12.py evt_400.fits --bands=c > v12_400.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/401     
# nohup ../eRASS_det_pipe_V12.py evt_401.fits --bands=c > v12_401.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/402     
# nohup ../eRASS_det_pipe_V12.py evt_402.fits --bands=c > v12_402.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/403     
# nohup ../eRASS_det_pipe_V12.py evt_403.fits --bands=c > v12_403.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/404     
# nohup ../eRASS_det_pipe_V12.py evt_404.fits --bands=c > v12_404.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/405     
# nohup ../eRASS_det_pipe_V12.py evt_405.fits --bands=c > v12_405.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/406     
# nohup ../eRASS_det_pipe_V12.py evt_406.fits --bands=c > v12_406.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/407     
# nohup ../eRASS_det_pipe_V12.py evt_407.fits --bands=c > v12_407.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/408     
# nohup ../eRASS_det_pipe_V12.py evt_408.fits --bands=c > v12_408.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/409     
# nohup ../eRASS_det_pipe_V12.py evt_409.fits --bands=c > v12_409.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/410     
# nohup ../eRASS_det_pipe_V12.py evt_410.fits --bands=c > v12_410.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/411     
# nohup ../eRASS_det_pipe_V12.py evt_411.fits --bands=c > v12_411.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/412     
# nohup ../eRASS_det_pipe_V12.py evt_412.fits --bands=c > v12_412.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/413     
# nohup ../eRASS_det_pipe_V12.py evt_413.fits --bands=c > v12_413.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/414     
# nohup ../eRASS_det_pipe_V12.py evt_414.fits --bands=c > v12_414.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/415     
# nohup ../eRASS_det_pipe_V12.py evt_415.fits --bands=c > v12_415.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/416     
# nohup ../eRASS_det_pipe_V12.py evt_416.fits --bands=c > v12_416.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/417     
# nohup ../eRASS_det_pipe_V12.py evt_417.fits --bands=c > v12_417.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/418     
# nohup ../eRASS_det_pipe_V12.py evt_418.fits --bands=c > v12_418.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/419     
# nohup ../eRASS_det_pipe_V12.py evt_419.fits --bands=c > v12_419.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/420     
# nohup ../eRASS_det_pipe_V12.py evt_420.fits --bands=c > v12_420.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/421     
# nohup ../eRASS_det_pipe_V12.py evt_421.fits --bands=c > v12_421.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/422     
# nohup ../eRASS_det_pipe_V12.py evt_422.fits --bands=c > v12_422.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/423     
# nohup ../eRASS_det_pipe_V12.py evt_423.fits --bands=c > v12_423.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/424     
# nohup ../eRASS_det_pipe_V12.py evt_424.fits --bands=c > v12_424.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/425     
# nohup ../eRASS_det_pipe_V12.py evt_425.fits --bands=c > v12_425.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/426     
# nohup ../eRASS_det_pipe_V12.py evt_426.fits --bands=c > v12_426.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/427     
# nohup ../eRASS_det_pipe_V12.py evt_427.fits --bands=c > v12_427.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/428     
# nohup ../eRASS_det_pipe_V12.py evt_428.fits --bands=c > v12_428.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/429     
# nohup ../eRASS_det_pipe_V12.py evt_429.fits --bands=c > v12_429.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/430     
# nohup ../eRASS_det_pipe_V12.py evt_430.fits --bands=c > v12_430.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/431     
# nohup ../eRASS_det_pipe_V12.py evt_431.fits --bands=c > v12_431.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/432     
# nohup ../eRASS_det_pipe_V12.py evt_432.fits --bands=c > v12_432.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/433     
# nohup ../eRASS_det_pipe_V12.py evt_433.fits --bands=c > v12_433.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/434     
# nohup ../eRASS_det_pipe_V12.py evt_434.fits --bands=c > v12_434.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/435     
# nohup ../eRASS_det_pipe_V12.py evt_435.fits --bands=c > v12_435.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/436     
# nohup ../eRASS_det_pipe_V12.py evt_436.fits --bands=c > v12_436.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/437     
# nohup ../eRASS_det_pipe_V12.py evt_437.fits --bands=c > v12_437.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/438     
# nohup ../eRASS_det_pipe_V12.py evt_438.fits --bands=c > v12_438.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/439     
# nohup ../eRASS_det_pipe_V12.py evt_439.fits --bands=c > v12_439.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/440     
# nohup ../eRASS_det_pipe_V12.py evt_440.fits --bands=c > v12_440.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/441     
# nohup ../eRASS_det_pipe_V12.py evt_441.fits --bands=c > v12_441.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/442     
# nohup ../eRASS_det_pipe_V12.py evt_442.fits --bands=c > v12_442.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/443     
# nohup ../eRASS_det_pipe_V12.py evt_443.fits --bands=c > v12_443.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/444     
# nohup ../eRASS_det_pipe_V12.py evt_444.fits --bands=c > v12_444.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/445     
# nohup ../eRASS_det_pipe_V12.py evt_445.fits --bands=c > v12_445.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/446     
# nohup ../eRASS_det_pipe_V12.py evt_446.fits --bands=c > v12_446.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/447     
# nohup ../eRASS_det_pipe_V12.py evt_447.fits --bands=c > v12_447.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/448     
# nohup ../eRASS_det_pipe_V12.py evt_448.fits --bands=c > v12_448.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/449     
# nohup ../eRASS_det_pipe_V12.py evt_449.fits --bands=c > v12_449.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/450     
# nohup ../eRASS_det_pipe_V12.py evt_450.fits --bands=c > v12_450.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/451     
# nohup ../eRASS_det_pipe_V12.py evt_451.fits --bands=c > v12_451.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/452     
# nohup ../eRASS_det_pipe_V12.py evt_452.fits --bands=c > v12_452.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/453     
# nohup ../eRASS_det_pipe_V12.py evt_453.fits --bands=c > v12_453.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/454     
# nohup ../eRASS_det_pipe_V12.py evt_454.fits --bands=c > v12_454.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/455     
# nohup ../eRASS_det_pipe_V12.py evt_455.fits --bands=c > v12_455.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/456     
# nohup ../eRASS_det_pipe_V12.py evt_456.fits --bands=c > v12_456.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/457     
# nohup ../eRASS_det_pipe_V12.py evt_457.fits --bands=c > v12_457.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/458     
# nohup ../eRASS_det_pipe_V12.py evt_458.fits --bands=c > v12_458.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/459     
# nohup ../eRASS_det_pipe_V12.py evt_459.fits --bands=c > v12_459.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/460     
# nohup ../eRASS_det_pipe_V12.py evt_460.fits --bands=c > v12_460.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/461     
# nohup ../eRASS_det_pipe_V12.py evt_461.fits --bands=c > v12_461.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/462     
# nohup ../eRASS_det_pipe_V12.py evt_462.fits --bands=c > v12_462.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/463     
# nohup ../eRASS_det_pipe_V12.py evt_463.fits --bands=c > v12_463.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/464     
# nohup ../eRASS_det_pipe_V12.py evt_464.fits --bands=c > v12_464.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/465     
# nohup ../eRASS_det_pipe_V12.py evt_465.fits --bands=c > v12_465.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/466     
# nohup ../eRASS_det_pipe_V12.py evt_466.fits --bands=c > v12_466.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/467     
# nohup ../eRASS_det_pipe_V12.py evt_467.fits --bands=c > v12_467.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/468     
# nohup ../eRASS_det_pipe_V12.py evt_468.fits --bands=c > v12_468.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/469     
# nohup ../eRASS_det_pipe_V12.py evt_469.fits --bands=c > v12_469.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/470     
# nohup ../eRASS_det_pipe_V12.py evt_470.fits --bands=c > v12_470.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/471     
# nohup ../eRASS_det_pipe_V12.py evt_471.fits --bands=c > v12_471.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/472     
# nohup ../eRASS_det_pipe_V12.py evt_472.fits --bands=c > v12_472.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/473     
# nohup ../eRASS_det_pipe_V12.py evt_473.fits --bands=c > v12_473.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/474     
# nohup ../eRASS_det_pipe_V12.py evt_474.fits --bands=c > v12_474.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/475     
# nohup ../eRASS_det_pipe_V12.py evt_475.fits --bands=c > v12_475.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/476     
# nohup ../eRASS_det_pipe_V12.py evt_476.fits --bands=c > v12_476.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/477     
# nohup ../eRASS_det_pipe_V12.py evt_477.fits --bands=c > v12_477.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/478     
# nohup ../eRASS_det_pipe_V12.py evt_478.fits --bands=c > v12_478.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/479     
# nohup ../eRASS_det_pipe_V12.py evt_479.fits --bands=c > v12_479.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/480     
# nohup ../eRASS_det_pipe_V12.py evt_480.fits --bands=c > v12_480.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/481     
nohup ../eRASS_det_pipe_V12.py evt_481.fits --bands=c > v12_481.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/482     
nohup ../eRASS_det_pipe_V12.py evt_482.fits --bands=c > v12_482.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/483     
nohup ../eRASS_det_pipe_V12.py evt_483.fits --bands=c > v12_483.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/484     
nohup ../eRASS_det_pipe_V12.py evt_484.fits --bands=c > v12_484.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/485     
nohup ../eRASS_det_pipe_V12.py evt_485.fits --bands=c > v12_485.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/486     
nohup ../eRASS_det_pipe_V12.py evt_486.fits --bands=c > v12_486.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/487     
nohup ../eRASS_det_pipe_V12.py evt_487.fits --bands=c > v12_487.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/488     
nohup ../eRASS_det_pipe_V12.py evt_488.fits --bands=c > v12_488.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/489     
nohup ../eRASS_det_pipe_V12.py evt_489.fits --bands=c > v12_489.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/490     
nohup ../eRASS_det_pipe_V12.py evt_490.fits --bands=c > v12_490.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/491     
# nohup ../eRASS_det_pipe_V12.py evt_491.fits --bands=c > v12_491.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/492     
# nohup ../eRASS_det_pipe_V12.py evt_492.fits --bands=c > v12_492.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/493     
# nohup ../eRASS_det_pipe_V12.py evt_493.fits --bands=c > v12_493.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/494     
# nohup ../eRASS_det_pipe_V12.py evt_494.fits --bands=c > v12_494.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/495     
# nohup ../eRASS_det_pipe_V12.py evt_495.fits --bands=c > v12_495.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/496     
# nohup ../eRASS_det_pipe_V12.py evt_496.fits --bands=c > v12_496.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/497     
# nohup ../eRASS_det_pipe_V12.py evt_497.fits --bands=c > v12_497.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/498     
# nohup ../eRASS_det_pipe_V12.py evt_498.fits --bands=c > v12_498.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/499     
# nohup ../eRASS_det_pipe_V12.py evt_499.fits --bands=c > v12_499.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/500     
# nohup ../eRASS_det_pipe_V12.py evt_500.fits --bands=c > v12_500.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/501     
# nohup ../eRASS_det_pipe_V12.py evt_501.fits --bands=c > v12_501.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/502     
# nohup ../eRASS_det_pipe_V12.py evt_502.fits --bands=c > v12_502.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/503     
# nohup ../eRASS_det_pipe_V12.py evt_503.fits --bands=c > v12_503.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/504     
# nohup ../eRASS_det_pipe_V12.py evt_504.fits --bands=c > v12_504.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/505     
# nohup ../eRASS_det_pipe_V12.py evt_505.fits --bands=c > v12_505.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/506     
# nohup ../eRASS_det_pipe_V12.py evt_506.fits --bands=c > v12_506.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/507     
# nohup ../eRASS_det_pipe_V12.py evt_507.fits --bands=c > v12_507.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/508     
# nohup ../eRASS_det_pipe_V12.py evt_508.fits --bands=c > v12_508.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/509     
# nohup ../eRASS_det_pipe_V12.py evt_509.fits --bands=c > v12_509.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/510     
nohup ../eRASS_det_pipe_V12.py evt_510.fits --bands=c > v12_510.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/511     
# nohup ../eRASS_det_pipe_V12.py evt_511.fits --bands=c > v12_511.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/512     
# nohup ../eRASS_det_pipe_V12.py evt_512.fits --bands=c > v12_512.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/513     
# nohup ../eRASS_det_pipe_V12.py evt_513.fits --bands=c > v12_513.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/514     
# nohup ../eRASS_det_pipe_V12.py evt_514.fits --bands=c > v12_514.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/515     
# nohup ../eRASS_det_pipe_V12.py evt_515.fits --bands=c > v12_515.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/516     
nohup ../eRASS_det_pipe_V12.py evt_516.fits --bands=c > v12_516.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/517     
nohup ../eRASS_det_pipe_V12.py evt_517.fits --bands=c > v12_517.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/518     
nohup ../eRASS_det_pipe_V12.py evt_518.fits --bands=c > v12_518.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/519     
nohup ../eRASS_det_pipe_V12.py evt_519.fits --bands=c > v12_519.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/520     
# nohup ../eRASS_det_pipe_V12.py evt_520.fits --bands=c > v12_520.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/521     
# nohup ../eRASS_det_pipe_V12.py evt_521.fits --bands=c > v12_521.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/522     
# nohup ../eRASS_det_pipe_V12.py evt_522.fits --bands=c > v12_522.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/523     
# nohup ../eRASS_det_pipe_V12.py evt_523.fits --bands=c > v12_523.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/524     
# nohup ../eRASS_det_pipe_V12.py evt_524.fits --bands=c > v12_524.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/525     
# nohup ../eRASS_det_pipe_V12.py evt_525.fits --bands=c > v12_525.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/526     
# nohup ../eRASS_det_pipe_V12.py evt_526.fits --bands=c > v12_526.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/527     
# nohup ../eRASS_det_pipe_V12.py evt_527.fits --bands=c > v12_527.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/528     
nohup ../eRASS_det_pipe_V12.py evt_528.fits --bands=c > v12_528.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/529     
nohup ../eRASS_det_pipe_V12.py evt_529.fits --bands=c > v12_529.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/530     
nohup ../eRASS_det_pipe_V12.py evt_530.fits --bands=c > v12_530.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/531     
nohup ../eRASS_det_pipe_V12.py evt_531.fits --bands=c > v12_531.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/532     
nohup ../eRASS_det_pipe_V12.py evt_532.fits --bands=c > v12_532.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/533     
# nohup ../eRASS_det_pipe_V12.py evt_533.fits --bands=c > v12_533.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/534     
# nohup ../eRASS_det_pipe_V12.py evt_534.fits --bands=c > v12_534.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/535     
# nohup ../eRASS_det_pipe_V12.py evt_535.fits --bands=c > v12_535.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/536     
# nohup ../eRASS_det_pipe_V12.py evt_536.fits --bands=c > v12_536.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/537     
# nohup ../eRASS_det_pipe_V12.py evt_537.fits --bands=c > v12_537.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/538     
# nohup ../eRASS_det_pipe_V12.py evt_538.fits --bands=c > v12_538.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/539     
# nohup ../eRASS_det_pipe_V12.py evt_539.fits --bands=c > v12_539.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/540     
# nohup ../eRASS_det_pipe_V12.py evt_540.fits --bands=c > v12_540.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/541     
# nohup ../eRASS_det_pipe_V12.py evt_541.fits --bands=c > v12_541.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/542     
# nohup ../eRASS_det_pipe_V12.py evt_542.fits --bands=c > v12_542.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/543     
# nohup ../eRASS_det_pipe_V12.py evt_543.fits --bands=c > v12_543.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/544     
# nohup ../eRASS_det_pipe_V12.py evt_544.fits --bands=c > v12_544.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/545     
# nohup ../eRASS_det_pipe_V12.py evt_545.fits --bands=c > v12_545.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/546     
# nohup ../eRASS_det_pipe_V12.py evt_546.fits --bands=c > v12_546.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/547     
# nohup ../eRASS_det_pipe_V12.py evt_547.fits --bands=c > v12_547.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/548     
# nohup ../eRASS_det_pipe_V12.py evt_548.fits --bands=c > v12_548.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/549     
# nohup ../eRASS_det_pipe_V12.py evt_549.fits --bands=c > v12_549.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/550     
# nohup ../eRASS_det_pipe_V12.py evt_550.fits --bands=c > v12_550.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/551     
# nohup ../eRASS_det_pipe_V12.py evt_551.fits --bands=c > v12_551.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/552     
# nohup ../eRASS_det_pipe_V12.py evt_552.fits --bands=c > v12_552.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/553     
# nohup ../eRASS_det_pipe_V12.py evt_553.fits --bands=c > v12_553.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/554     
# nohup ../eRASS_det_pipe_V12.py evt_554.fits --bands=c > v12_554.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/555     
# nohup ../eRASS_det_pipe_V12.py evt_555.fits --bands=c > v12_555.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/556     
# nohup ../eRASS_det_pipe_V12.py evt_556.fits --bands=c > v12_556.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/557     
# nohup ../eRASS_det_pipe_V12.py evt_557.fits --bands=c > v12_557.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/558     
# nohup ../eRASS_det_pipe_V12.py evt_558.fits --bands=c > v12_558.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/559     
# nohup ../eRASS_det_pipe_V12.py evt_559.fits --bands=c > v12_559.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/560     
# nohup ../eRASS_det_pipe_V12.py evt_560.fits --bands=c > v12_560.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/561     
# nohup ../eRASS_det_pipe_V12.py evt_561.fits --bands=c > v12_561.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/562     
# nohup ../eRASS_det_pipe_V12.py evt_562.fits --bands=c > v12_562.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/563     
# nohup ../eRASS_det_pipe_V12.py evt_563.fits --bands=c > v12_563.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/564     
# nohup ../eRASS_det_pipe_V12.py evt_564.fits --bands=c > v12_564.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/565     
# nohup ../eRASS_det_pipe_V12.py evt_565.fits --bands=c > v12_565.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/566     
# nohup ../eRASS_det_pipe_V12.py evt_566.fits --bands=c > v12_566.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/567     
# nohup ../eRASS_det_pipe_V12.py evt_567.fits --bands=c > v12_567.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/568     
# nohup ../eRASS_det_pipe_V12.py evt_568.fits --bands=c > v12_568.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/569     
# nohup ../eRASS_det_pipe_V12.py evt_569.fits --bands=c > v12_569.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/570     
# nohup ../eRASS_det_pipe_V12.py evt_570.fits --bands=c > v12_570.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/571     
# nohup ../eRASS_det_pipe_V12.py evt_571.fits --bands=c > v12_571.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/572     
# nohup ../eRASS_det_pipe_V12.py evt_572.fits --bands=c > v12_572.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/573     
# nohup ../eRASS_det_pipe_V12.py evt_573.fits --bands=c > v12_573.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/574     
# nohup ../eRASS_det_pipe_V12.py evt_574.fits --bands=c > v12_574.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/575     
# nohup ../eRASS_det_pipe_V12.py evt_575.fits --bands=c > v12_575.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/576     
# nohup ../eRASS_det_pipe_V12.py evt_576.fits --bands=c > v12_576.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/577     
# nohup ../eRASS_det_pipe_V12.py evt_577.fits --bands=c > v12_577.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/578     
# nohup ../eRASS_det_pipe_V12.py evt_578.fits --bands=c > v12_578.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/579     
# nohup ../eRASS_det_pipe_V12.py evt_579.fits --bands=c > v12_579.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/580     
# nohup ../eRASS_det_pipe_V12.py evt_580.fits --bands=c > v12_580.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/581     
# nohup ../eRASS_det_pipe_V12.py evt_581.fits --bands=c > v12_581.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/582     
# nohup ../eRASS_det_pipe_V12.py evt_582.fits --bands=c > v12_582.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/583     
# nohup ../eRASS_det_pipe_V12.py evt_583.fits --bands=c > v12_583.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/584     
# nohup ../eRASS_det_pipe_V12.py evt_584.fits --bands=c > v12_584.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/585     
# nohup ../eRASS_det_pipe_V12.py evt_585.fits --bands=c > v12_585.log & 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/586     
nohup ../eRASS_det_pipe_V12.py evt_586.fits --bands=c > v12_586.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/587     
# nohup ../eRASS_det_pipe_V12.py evt_587.fits --bands=c > v12_587.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/588     
# nohup ../eRASS_det_pipe_V12.py evt_588.fits --bands=c > v12_588.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/589     
# nohup ../eRASS_det_pipe_V12.py evt_589.fits --bands=c > v12_589.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/590     
# nohup ../eRASS_det_pipe_V12.py evt_590.fits --bands=c > v12_590.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/591     
# nohup ../eRASS_det_pipe_V12.py evt_591.fits --bands=c > v12_591.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/592     
# nohup ../eRASS_det_pipe_V12.py evt_592.fits --bands=c > v12_592.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/593     
# nohup ../eRASS_det_pipe_V12.py evt_593.fits --bands=c > v12_593.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/594     
# nohup ../eRASS_det_pipe_V12.py evt_594.fits --bands=c > v12_594.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/595     
# nohup ../eRASS_det_pipe_V12.py evt_595.fits --bands=c > v12_595.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/596     
# nohup ../eRASS_det_pipe_V12.py evt_596.fits --bands=c > v12_596.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/597     
# nohup ../eRASS_det_pipe_V12.py evt_597.fits --bands=c > v12_597.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/598     
# nohup ../eRASS_det_pipe_V12.py evt_598.fits --bands=c > v12_598.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/599     
# nohup ../eRASS_det_pipe_V12.py evt_599.fits --bands=c > v12_599.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/600     
# nohup ../eRASS_det_pipe_V12.py evt_600.fits --bands=c > v12_600.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/601     
# nohup ../eRASS_det_pipe_V12.py evt_601.fits --bands=c > v12_601.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/602     
# nohup ../eRASS_det_pipe_V12.py evt_602.fits --bands=c > v12_602.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/603     
# nohup ../eRASS_det_pipe_V12.py evt_603.fits --bands=c > v12_603.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/604     
# nohup ../eRASS_det_pipe_V12.py evt_604.fits --bands=c > v12_604.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/605     
# nohup ../eRASS_det_pipe_V12.py evt_605.fits --bands=c > v12_605.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/606     
# nohup ../eRASS_det_pipe_V12.py evt_606.fits --bands=c > v12_606.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/607     
# nohup ../eRASS_det_pipe_V12.py evt_607.fits --bands=c > v12_607.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/608     
# nohup ../eRASS_det_pipe_V12.py evt_608.fits --bands=c > v12_608.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/609     
# nohup ../eRASS_det_pipe_V12.py evt_609.fits --bands=c > v12_609.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/610     
# nohup ../eRASS_det_pipe_V12.py evt_610.fits --bands=c > v12_610.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/611     
# nohup ../eRASS_det_pipe_V12.py evt_611.fits --bands=c > v12_611.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/612     
# nohup ../eRASS_det_pipe_V12.py evt_612.fits --bands=c > v12_612.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/613     
# nohup ../eRASS_det_pipe_V12.py evt_613.fits --bands=c > v12_613.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/614     
# nohup ../eRASS_det_pipe_V12.py evt_614.fits --bands=c > v12_614.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/615     
# nohup ../eRASS_det_pipe_V12.py evt_615.fits --bands=c > v12_615.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/616     
# nohup ../eRASS_det_pipe_V12.py evt_616.fits --bands=c > v12_616.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/617     
# nohup ../eRASS_det_pipe_V12.py evt_617.fits --bands=c > v12_617.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/618     
# nohup ../eRASS_det_pipe_V12.py evt_618.fits --bands=c > v12_618.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/619     
# nohup ../eRASS_det_pipe_V12.py evt_619.fits --bands=c > v12_619.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/620     
# nohup ../eRASS_det_pipe_V12.py evt_620.fits --bands=c > v12_620.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/621     
# nohup ../eRASS_det_pipe_V12.py evt_621.fits --bands=c > v12_621.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/622     
# nohup ../eRASS_det_pipe_V12.py evt_622.fits --bands=c > v12_622.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/623     
# nohup ../eRASS_det_pipe_V12.py evt_623.fits --bands=c > v12_623.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/624     
# nohup ../eRASS_det_pipe_V12.py evt_624.fits --bands=c > v12_624.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/625     
# nohup ../eRASS_det_pipe_V12.py evt_625.fits --bands=c > v12_625.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/626     
# nohup ../eRASS_det_pipe_V12.py evt_626.fits --bands=c > v12_626.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/627     
# nohup ../eRASS_det_pipe_V12.py evt_627.fits --bands=c > v12_627.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/628     
# nohup ../eRASS_det_pipe_V12.py evt_628.fits --bands=c > v12_628.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/629     
# nohup ../eRASS_det_pipe_V12.py evt_629.fits --bands=c > v12_629.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/630     
# nohup ../eRASS_det_pipe_V12.py evt_630.fits --bands=c > v12_630.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/631     
# nohup ../eRASS_det_pipe_V12.py evt_631.fits --bands=c > v12_631.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/632     
# nohup ../eRASS_det_pipe_V12.py evt_632.fits --bands=c > v12_632.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/633     
# nohup ../eRASS_det_pipe_V12.py evt_633.fits --bands=c > v12_633.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/634     
# nohup ../eRASS_det_pipe_V12.py evt_634.fits --bands=c > v12_634.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/635     
# nohup ../eRASS_det_pipe_V12.py evt_635.fits --bands=c > v12_635.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/636     
# nohup ../eRASS_det_pipe_V12.py evt_636.fits --bands=c > v12_636.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/637     
# nohup ../eRASS_det_pipe_V12.py evt_637.fits --bands=c > v12_637.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/638     
# nohup ../eRASS_det_pipe_V12.py evt_638.fits --bands=c > v12_638.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/639     
# nohup ../eRASS_det_pipe_V12.py evt_639.fits --bands=c > v12_639.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/640     
# nohup ../eRASS_det_pipe_V12.py evt_640.fits --bands=c > v12_640.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/641     
# nohup ../eRASS_det_pipe_V12.py evt_641.fits --bands=c > v12_641.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/642     
# nohup ../eRASS_det_pipe_V12.py evt_642.fits --bands=c > v12_642.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/643     
# nohup ../eRASS_det_pipe_V12.py evt_643.fits --bands=c > v12_643.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/644     
# nohup ../eRASS_det_pipe_V12.py evt_644.fits --bands=c > v12_644.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/645     
# nohup ../eRASS_det_pipe_V12.py evt_645.fits --bands=c > v12_645.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/646     
# nohup ../eRASS_det_pipe_V12.py evt_646.fits --bands=c > v12_646.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/647     
# nohup ../eRASS_det_pipe_V12.py evt_647.fits --bands=c > v12_647.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/648     
# nohup ../eRASS_det_pipe_V12.py evt_648.fits --bands=c > v12_648.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/649     
# nohup ../eRASS_det_pipe_V12.py evt_649.fits --bands=c > v12_649.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/650     
# nohup ../eRASS_det_pipe_V12.py evt_650.fits --bands=c > v12_650.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/651     
# nohup ../eRASS_det_pipe_V12.py evt_651.fits --bands=c > v12_651.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/652     
# nohup ../eRASS_det_pipe_V12.py evt_652.fits --bands=c > v12_652.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/653     
# nohup ../eRASS_det_pipe_V12.py evt_653.fits --bands=c > v12_653.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/654     
# nohup ../eRASS_det_pipe_V12.py evt_654.fits --bands=c > v12_654.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/655     
# nohup ../eRASS_det_pipe_V12.py evt_655.fits --bands=c > v12_655.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/656     
# nohup ../eRASS_det_pipe_V12.py evt_656.fits --bands=c > v12_656.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/657     
# nohup ../eRASS_det_pipe_V12.py evt_657.fits --bands=c > v12_657.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/658     
# nohup ../eRASS_det_pipe_V12.py evt_658.fits --bands=c > v12_658.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/659     
# nohup ../eRASS_det_pipe_V12.py evt_659.fits --bands=c > v12_659.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/660     
# nohup ../eRASS_det_pipe_V12.py evt_660.fits --bands=c > v12_660.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/661     
# nohup ../eRASS_det_pipe_V12.py evt_661.fits --bands=c > v12_661.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/662     
# nohup ../eRASS_det_pipe_V12.py evt_662.fits --bands=c > v12_662.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/663     
# nohup ../eRASS_det_pipe_V12.py evt_663.fits --bands=c > v12_663.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/664     
# nohup ../eRASS_det_pipe_V12.py evt_664.fits --bands=c > v12_664.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/665     
# nohup ../eRASS_det_pipe_V12.py evt_665.fits --bands=c > v12_665.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/666     
# nohup ../eRASS_det_pipe_V12.py evt_666.fits --bands=c > v12_666.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/667     
# nohup ../eRASS_det_pipe_V12.py evt_667.fits --bands=c > v12_667.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/668     
# nohup ../eRASS_det_pipe_V12.py evt_668.fits --bands=c > v12_668.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/669     
# nohup ../eRASS_det_pipe_V12.py evt_669.fits --bands=c > v12_669.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/670     
# nohup ../eRASS_det_pipe_V12.py evt_670.fits --bands=c > v12_670.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/671     
# nohup ../eRASS_det_pipe_V12.py evt_671.fits --bands=c > v12_671.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/672     
# nohup ../eRASS_det_pipe_V12.py evt_672.fits --bands=c > v12_672.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/673     
# nohup ../eRASS_det_pipe_V12.py evt_673.fits --bands=c > v12_673.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/674     
# nohup ../eRASS_det_pipe_V12.py evt_674.fits --bands=c > v12_674.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/675     
# nohup ../eRASS_det_pipe_V12.py evt_675.fits --bands=c > v12_675.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/676     
# nohup ../eRASS_det_pipe_V12.py evt_676.fits --bands=c > v12_676.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/677     
# nohup ../eRASS_det_pipe_V12.py evt_677.fits --bands=c > v12_677.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/678     
# nohup ../eRASS_det_pipe_V12.py evt_678.fits --bands=c > v12_678.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/679     
# nohup ../eRASS_det_pipe_V12.py evt_679.fits --bands=c > v12_679.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/680     
# nohup ../eRASS_det_pipe_V12.py evt_680.fits --bands=c > v12_680.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/681     
# nohup ../eRASS_det_pipe_V12.py evt_681.fits --bands=c > v12_681.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/682     
# nohup ../eRASS_det_pipe_V12.py evt_682.fits --bands=c > v12_682.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/683     
# nohup ../eRASS_det_pipe_V12.py evt_683.fits --bands=c > v12_683.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/684     
# nohup ../eRASS_det_pipe_V12.py evt_684.fits --bands=c > v12_684.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/685     
# nohup ../eRASS_det_pipe_V12.py evt_685.fits --bands=c > v12_685.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/686     
# nohup ../eRASS_det_pipe_V12.py evt_686.fits --bands=c > v12_686.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/687     
# nohup ../eRASS_det_pipe_V12.py evt_687.fits --bands=c > v12_687.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/688     
# nohup ../eRASS_det_pipe_V12.py evt_688.fits --bands=c > v12_688.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/689     
# nohup ../eRASS_det_pipe_V12.py evt_689.fits --bands=c > v12_689.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/690     
# nohup ../eRASS_det_pipe_V12.py evt_690.fits --bands=c > v12_690.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/691     
# nohup ../eRASS_det_pipe_V12.py evt_691.fits --bands=c > v12_691.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/692     
# nohup ../eRASS_det_pipe_V12.py evt_692.fits --bands=c > v12_692.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/693     
# nohup ../eRASS_det_pipe_V12.py evt_693.fits --bands=c > v12_693.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/694     
# nohup ../eRASS_det_pipe_V12.py evt_694.fits --bands=c > v12_694.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/695     
# nohup ../eRASS_det_pipe_V12.py evt_695.fits --bands=c > v12_695.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/696     
# nohup ../eRASS_det_pipe_V12.py evt_696.fits --bands=c > v12_696.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/697     
# nohup ../eRASS_det_pipe_V12.py evt_697.fits --bands=c > v12_697.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/698     
# nohup ../eRASS_det_pipe_V12.py evt_698.fits --bands=c > v12_698.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/699     
# nohup ../eRASS_det_pipe_V12.py evt_699.fits --bands=c > v12_699.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/700     
# nohup ../eRASS_det_pipe_V12.py evt_700.fits --bands=c > v12_700.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/701     
# nohup ../eRASS_det_pipe_V12.py evt_701.fits --bands=c > v12_701.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/702     
# nohup ../eRASS_det_pipe_V12.py evt_702.fits --bands=c > v12_702.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/703     
# nohup ../eRASS_det_pipe_V12.py evt_703.fits --bands=c > v12_703.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/704     
# nohup ../eRASS_det_pipe_V12.py evt_704.fits --bands=c > v12_704.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/705     
# nohup ../eRASS_det_pipe_V12.py evt_705.fits --bands=c > v12_705.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/706     
# nohup ../eRASS_det_pipe_V12.py evt_706.fits --bands=c > v12_706.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/707     
# nohup ../eRASS_det_pipe_V12.py evt_707.fits --bands=c > v12_707.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/708     
# nohup ../eRASS_det_pipe_V12.py evt_708.fits --bands=c > v12_708.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/709     
# nohup ../eRASS_det_pipe_V12.py evt_709.fits --bands=c > v12_709.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/710     
# nohup ../eRASS_det_pipe_V12.py evt_710.fits --bands=c > v12_710.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/711     
# nohup ../eRASS_det_pipe_V12.py evt_711.fits --bands=c > v12_711.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/712     
# nohup ../eRASS_det_pipe_V12.py evt_712.fits --bands=c > v12_712.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/713     
# nohup ../eRASS_det_pipe_V12.py evt_713.fits --bands=c > v12_713.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/714     
# nohup ../eRASS_det_pipe_V12.py evt_714.fits --bands=c > v12_714.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/715     
# nohup ../eRASS_det_pipe_V12.py evt_715.fits --bands=c > v12_715.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/716     
# nohup ../eRASS_det_pipe_V12.py evt_716.fits --bands=c > v12_716.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/717     
# nohup ../eRASS_det_pipe_V12.py evt_717.fits --bands=c > v12_717.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/718     
# nohup ../eRASS_det_pipe_V12.py evt_718.fits --bands=c > v12_718.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/719     
# nohup ../eRASS_det_pipe_V12.py evt_719.fits --bands=c > v12_719.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/720     
# nohup ../eRASS_det_pipe_V12.py evt_720.fits --bands=c > v12_720.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/721     
# nohup ../eRASS_det_pipe_V12.py evt_721.fits --bands=c > v12_721.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/722     
# nohup ../eRASS_det_pipe_V12.py evt_722.fits --bands=c > v12_722.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/723     
# nohup ../eRASS_det_pipe_V12.py evt_723.fits --bands=c > v12_723.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/724     
# nohup ../eRASS_det_pipe_V12.py evt_724.fits --bands=c > v12_724.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/725     
# nohup ../eRASS_det_pipe_V12.py evt_725.fits --bands=c > v12_725.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/726     
# nohup ../eRASS_det_pipe_V12.py evt_726.fits --bands=c > v12_726.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/727     
# nohup ../eRASS_det_pipe_V12.py evt_727.fits --bands=c > v12_727.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/728     
# nohup ../eRASS_det_pipe_V12.py evt_728.fits --bands=c > v12_728.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/729     
# nohup ../eRASS_det_pipe_V12.py evt_729.fits --bands=c > v12_729.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/730     
# nohup ../eRASS_det_pipe_V12.py evt_730.fits --bands=c > v12_730.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/731     
# nohup ../eRASS_det_pipe_V12.py evt_731.fits --bands=c > v12_731.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/732     
# nohup ../eRASS_det_pipe_V12.py evt_732.fits --bands=c > v12_732.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/733     
# nohup ../eRASS_det_pipe_V12.py evt_733.fits --bands=c > v12_733.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/734     
# nohup ../eRASS_det_pipe_V12.py evt_734.fits --bands=c > v12_734.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/735     
# nohup ../eRASS_det_pipe_V12.py evt_735.fits --bands=c > v12_735.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/736     
# nohup ../eRASS_det_pipe_V12.py evt_736.fits --bands=c > v12_736.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/737     
# nohup ../eRASS_det_pipe_V12.py evt_737.fits --bands=c > v12_737.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/738     
# nohup ../eRASS_det_pipe_V12.py evt_738.fits --bands=c > v12_738.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/739     
# nohup ../eRASS_det_pipe_V12.py evt_739.fits --bands=c > v12_739.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/740     
# nohup ../eRASS_det_pipe_V12.py evt_740.fits --bands=c > v12_740.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/741     
# nohup ../eRASS_det_pipe_V12.py evt_741.fits --bands=c > v12_741.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/742     
# nohup ../eRASS_det_pipe_V12.py evt_742.fits --bands=c > v12_742.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/743     
# nohup ../eRASS_det_pipe_V12.py evt_743.fits --bands=c > v12_743.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/744     
# nohup ../eRASS_det_pipe_V12.py evt_744.fits --bands=c > v12_744.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/745     
# nohup ../eRASS_det_pipe_V12.py evt_745.fits --bands=c > v12_745.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/746     
# nohup ../eRASS_det_pipe_V12.py evt_746.fits --bands=c > v12_746.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/747     
# nohup ../eRASS_det_pipe_V12.py evt_747.fits --bands=c > v12_747.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/748     
# nohup ../eRASS_det_pipe_V12.py evt_748.fits --bands=c > v12_748.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/749     
# nohup ../eRASS_det_pipe_V12.py evt_749.fits --bands=c > v12_749.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/750     
# nohup ../eRASS_det_pipe_V12.py evt_750.fits --bands=c > v12_750.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/751     
# nohup ../eRASS_det_pipe_V12.py evt_751.fits --bands=c > v12_751.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/752     
# nohup ../eRASS_det_pipe_V12.py evt_752.fits --bands=c > v12_752.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/753     
# nohup ../eRASS_det_pipe_V12.py evt_753.fits --bands=c > v12_753.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/754     
# nohup ../eRASS_det_pipe_V12.py evt_754.fits --bands=c > v12_754.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/755     
# nohup ../eRASS_det_pipe_V12.py evt_755.fits --bands=c > v12_755.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/756     
# nohup ../eRASS_det_pipe_V12.py evt_756.fits --bands=c > v12_756.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/757     
# nohup ../eRASS_det_pipe_V12.py evt_757.fits --bands=c > v12_757.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/758     
# nohup ../eRASS_det_pipe_V12.py evt_758.fits --bands=c > v12_758.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/759     
# nohup ../eRASS_det_pipe_V12.py evt_759.fits --bands=c > v12_759.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/760     
# nohup ../eRASS_det_pipe_V12.py evt_760.fits --bands=c > v12_760.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/761     
# nohup ../eRASS_det_pipe_V12.py evt_761.fits --bands=c > v12_761.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/762     
# nohup ../eRASS_det_pipe_V12.py evt_762.fits --bands=c > v12_762.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/763     
# nohup ../eRASS_det_pipe_V12.py evt_763.fits --bands=c > v12_763.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/764     
# nohup ../eRASS_det_pipe_V12.py evt_764.fits --bands=c > v12_764.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/765     
# nohup ../eRASS_det_pipe_V12.py evt_765.fits --bands=c > v12_765.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/766     
# nohup ../eRASS_det_pipe_V12.py evt_766.fits --bands=c > v12_766.log & 
# cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/767     
# nohup ../eRASS_det_pipe_V12.py evt_767.fits --bands=c > v12_767.log & 

cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/201
../eRASS_det_pipe_V11.py evt_201.fits --bands=c > v12.log &
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/201
../eRASS_det_pipe_V11.py evt_201.fits --bands=c > v12.log &
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/201
../eRASS_det_pipe_V11.py evt_201.fits --bands=c > v12.log &
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/201
../eRASS_det_pipe_V11.py evt_201.fits --bands=c > v12.log &
 
pyCONDA
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/600
../eRASS_det_pipe_V11.py evt_600.fits --bands=c > v12.log &
 
pyCONDA
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/480
../eRASS_det_pipe_V11.py evt_480.fits --bands=c > v12.log &

pyCONDA
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/601
 ../eRASS_det_pipe_V11.py evt_601.fits --bands=c > v12.log &
 
pyCONDA
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
cd /data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn/600
 ../eRASS_det_pipe_V11.py evt_600.fits --bands=c > v12.log &

on ds52 

eSASSdevel              # EVTOOL bug
eSASSpipe_181017
eSASSusers_140526
eSASSusers_140905
eSASSusers_150701
eSASSusers_150901
eSASSusers_161028
eSASSusers_161130
eSASSusers_180416
eSASSusers_180416.saved
eSASSusers_190220-moved # EVTOOL bug
eSASSusers_190220       # seems to work !!
eSASSusers_190925-moved # EXPMAP bug
eSASSusers_190925       # seems to work !!
fits_090304

pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
nohup python esass.py 368 > esass_py_368.log & 
nohup python esass.py 400 > esass_py_400.log & 
nohup python esass.py 500 > esass_py_500.log & 

nohup python esass.py 000 > esass_py_000.log & 
nohup python esass.py 001 > esass_py_001.log & 
nohup python esass.py 002 > esass_py_002.log & 
nohup python esass.py 003 > esass_py_003.log & 
nohup python esass.py 004 > esass_py_004.log & 
nohup python esass.py 005 > esass_py_005.log & 
nohup python esass.py 006 > esass_py_006.log & 
nohup python esass.py 007 > esass_py_007.log & 
nohup python esass.py 008 > esass_py_008.log & 
nohup python esass.py 009 > esass_py_009.log & 
nohup python esass.py 010 > esass_py_010.log & 
nohup python esass.py 011 > esass_py_011.log & 
nohup python esass.py 012 > esass_py_012.log & 
nohup python esass.py 013 > esass_py_013.log & 
nohup python esass.py 014 > esass_py_014.log & 
nohup python esass.py 015 > esass_py_015.log & 
nohup python esass.py 016 > esass_py_016.log & 
nohup python esass.py 017 > esass_py_017.log & 
nohup python esass.py 018 > esass_py_018.log & 
nohup python esass.py 019 > esass_py_019.log & 
nohup python esass.py 020 > esass_py_020.log & 
nohup python esass.py 021 > esass_py_021.log & 
nohup python esass.py 022 > esass_py_022.log & 
nohup python esass.py 023 > esass_py_023.log & 
nohup python esass.py 024 > esass_py_024.log & 
nohup python esass.py 025 > esass_py_025.log & 
nohup python esass.py 026 > esass_py_026.log & 
nohup python esass.py 027 > esass_py_027.log & 
nohup python esass.py 028 > esass_py_028.log & 
nohup python esass.py 029 > esass_py_029.log & 
nohup python esass.py 030 > esass_py_030.log & 
nohup python esass.py 031 > esass_py_031.log & 
nohup python esass.py 032 > esass_py_032.log & 
nohup python esass.py 033 > esass_py_033.log & 
nohup python esass.py 034 > esass_py_034.log & 
nohup python esass.py 035 > esass_py_035.log & 
nohup python esass.py 036 > esass_py_036.log & 
nohup python esass.py 037 > esass_py_037.log & 
nohup python esass.py 038 > esass_py_038.log & 
nohup python esass.py 039 > esass_py_039.log & 
nohup python esass.py 040 > esass_py_040.log &
nohup python esass.py 041 > esass_py_041.log &
nohup python esass.py 042 > esass_py_042.log &
nohup python esass.py 043 > esass_py_043.log &
nohup python esass.py 044 > esass_py_044.log &
nohup python esass.py 045 > esass_py_045.log &
nohup python esass.py 046 > esass_py_046.log &
nohup python esass.py 047 > esass_py_047.log &
nohup python esass.py 048 > esass_py_048.log &
nohup python esass.py 049 > esass_py_049.log &
nohup python esass.py 050 > esass_py_050.log &
nohup python esass.py 051 > esass_py_051.log &
nohup python esass.py 052 > esass_py_052.log &
nohup python esass.py 053 > esass_py_053.log &
nohup python esass.py 054 > esass_py_054.log &
nohup python esass.py 055 > esass_py_055.log &
nohup python esass.py 056 > esass_py_056.log &
nohup python esass.py 057 > esass_py_057.log &
nohup python esass.py 058 > esass_py_058.log &
nohup python esass.py 059 > esass_py_059.log &
python esass.py 060 
python esass.py 061 
python esass.py 062 
python esass.py 063 
python esass.py 064 
python esass.py 065 
python esass.py 066 
python esass.py 067 
python esass.py 068 
python esass.py 069 
python esass.py 070 
python esass.py 071 
python esass.py 072 
python esass.py 073 
python esass.py 074 
python esass.py 075 
python esass.py 076 
python esass.py 077 
python esass.py 078 
python esass.py 079 
python esass.py 080 
python esass.py 081 
python esass.py 082 
python esass.py 083 
python esass.py 084 
python esass.py 085 
python esass.py 086 
python esass.py 087 
python esass.py 088 
python esass.py 089 
python esass.py 090 
python esass.py 091 
python esass.py 092 
python esass.py 093 
python esass.py 094 
python esass.py 095 
python esass.py 096 
python esass.py 097 
python esass.py 098 
nohup python esass.py 098 > esass_py_098.log & 
nohup python esass.py 099 > esass_py_099.log & 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 100 
python esass.py 101 
python esass.py 102 
python esass.py 103 
python esass.py 104 
python esass.py 105 
python esass.py 106 
python esass.py 107 
python esass.py 108 
python esass.py 109 
python esass.py 110 
python esass.py 111 
python esass.py 112 
python esass.py 113 
python esass.py 114 
python esass.py 115 
python esass.py 116 
python esass.py 117 
python esass.py 118 
python esass.py 119 
python esass.py 120 
python esass.py 121 
python esass.py 122 
python esass.py 123 
python esass.py 124 
python esass.py 125 
python esass.py 126 
python esass.py 127 
python esass.py 128 
python esass.py 129 
python esass.py 130 
python esass.py 131 
python esass.py 132 
python esass.py 133 
python esass.py 134 
python esass.py 135 
python esass.py 136 
python esass.py 137 
python esass.py 138 
python esass.py 139 
python esass.py 140 
python esass.py 141 
python esass.py 142 
python esass.py 143 
python esass.py 144 
python esass.py 145 
python esass.py 146 
python esass.py 147 
python esass.py 148 
python esass.py 149 
python esass.py 150 
python esass.py 151 
python esass.py 152 
python esass.py 153 
python esass.py 154 
python esass.py 155 
python esass.py 156 
python esass.py 157 
python esass.py 158 
python esass.py 159 
python esass.py 160 
python esass.py 161 
python esass.py 162 
python esass.py 163 
python esass.py 164 
python esass.py 165 
python esass.py 166 
python esass.py 167 
python esass.py 168 
python esass.py 169 
python esass.py 170 
python esass.py 171 
python esass.py 172 
python esass.py 173 
python esass.py 174 
python esass.py 175 
python esass.py 176 
python esass.py 177 
python esass.py 178 
python esass.py 179 
python esass.py 180 
python esass.py 181 
python esass.py 182 
python esass.py 183 
python esass.py 184 
python esass.py 185 
python esass.py 186 
python esass.py 187 
python esass.py 188 
python esass.py 189 
python esass.py 190 
python esass.py 191 
python esass.py 192 
python esass.py 193 
python esass.py 194 
python esass.py 195 
python esass.py 196 
python esass.py 197 
python esass.py 198 
python esass.py 199 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 200 
python esass.py 201 
python esass.py 202 
python esass.py 203 
python esass.py 204 
python esass.py 205 
python esass.py 206 
python esass.py 207 
python esass.py 208 
python esass.py 209 
python esass.py 210 
python esass.py 211 
python esass.py 212 
python esass.py 213 
python esass.py 214 
python esass.py 215 
python esass.py 216 
python esass.py 217 
python esass.py 218 
python esass.py 219 
python esass.py 220 
python esass.py 221 
python esass.py 222 
python esass.py 223 
python esass.py 224 
python esass.py 225 
python esass.py 226 
python esass.py 227 
python esass.py 228 
python esass.py 229 
python esass.py 230 
python esass.py 231 
python esass.py 232 
python esass.py 233 
python esass.py 234 
python esass.py 235 
python esass.py 236 
python esass.py 237 
python esass.py 238 
python esass.py 239 
python esass.py 240 
python esass.py 241 
python esass.py 242 
python esass.py 243 
python esass.py 244 
python esass.py 245 
python esass.py 246 
python esass.py 247 
python esass.py 248 
python esass.py 249 
python esass.py 250 
python esass.py 251 
python esass.py 252 
python esass.py 253 
python esass.py 254 
python esass.py 255 
python esass.py 256 
python esass.py 257 
python esass.py 258 
python esass.py 259 
python esass.py 260 
python esass.py 261 
python esass.py 262 
python esass.py 263 
python esass.py 264 
python esass.py 265 
python esass.py 266 
python esass.py 267 
python esass.py 268 
python esass.py 269 
python esass.py 270 
python esass.py 271 
python esass.py 272 
python esass.py 273 
python esass.py 274 
python esass.py 275 
python esass.py 276 
python esass.py 277 
python esass.py 278 
python esass.py 279 
python esass.py 280 
python esass.py 281 
python esass.py 282 
python esass.py 283 
python esass.py 284 
python esass.py 285 
python esass.py 286 
python esass.py 287 
python esass.py 288 
python esass.py 289 
python esass.py 290 
python esass.py 291 
python esass.py 292 
python esass.py 293 
python esass.py 294 
python esass.py 295 
python esass.py 296 
python esass.py 297 
python esass.py 298 
python esass.py 299 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 300 
python esass.py 301 
python esass.py 302 
python esass.py 303 
python esass.py 304 
python esass.py 305 
python esass.py 306 
python esass.py 307 
python esass.py 308 
python esass.py 309 
python esass.py 310 
python esass.py 311 
python esass.py 312 
python esass.py 313 
python esass.py 314 
python esass.py 315 
python esass.py 316 
python esass.py 317 
python esass.py 318 
python esass.py 319 
python esass.py 320 
python esass.py 321 
python esass.py 322 
python esass.py 323 
python esass.py 324 
python esass.py 325 
python esass.py 326 
python esass.py 327 
python esass.py 328 
python esass.py 329 
python esass.py 330 
python esass.py 331 
python esass.py 332 
python esass.py 333 
python esass.py 334 
python esass.py 335 
python esass.py 336 
python esass.py 337 
python esass.py 338 
python esass.py 339 
python esass.py 340 
python esass.py 341 
python esass.py 342 
python esass.py 343 
python esass.py 344 
python esass.py 345 
python esass.py 346 
python esass.py 347 
python esass.py 348 
python esass.py 349 
python esass.py 350 
python esass.py 351 
python esass.py 352 
python esass.py 353 
python esass.py 354 
python esass.py 355 
python esass.py 356 
python esass.py 357 
python esass.py 358 
python esass.py 359 
python esass.py 360 
python esass.py 361 
python esass.py 362 
python esass.py 363 
python esass.py 364 
python esass.py 365 
python esass.py 366 
python esass.py 367 
python esass.py 368 
python esass.py 369 
python esass.py 370 
python esass.py 371 
python esass.py 372 
python esass.py 373 
python esass.py 374 
python esass.py 375 
python esass.py 376 
python esass.py 377 
python esass.py 378 
python esass.py 379 
python esass.py 380 
python esass.py 381 
python esass.py 382 
python esass.py 383 
python esass.py 384 
python esass.py 385 
python esass.py 386 
python esass.py 387 
python esass.py 388 
python esass.py 389 
python esass.py 390 
python esass.py 391 
python esass.py 392 
python esass.py 393 
python esass.py 394 
python esass.py 395 
python esass.py 396 
python esass.py 397 
python esass.py 398 
python esass.py 399 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 400 
python esass.py 401 
python esass.py 402 
python esass.py 403 
python esass.py 404 
python esass.py 405 
python esass.py 406 
python esass.py 407 
python esass.py 408 
python esass.py 409 
python esass.py 410 
python esass.py 411 
python esass.py 412 
python esass.py 413 
python esass.py 414 
python esass.py 415 
python esass.py 416 
python esass.py 417 
python esass.py 418 
python esass.py 419 
python esass.py 420 
python esass.py 421 
python esass.py 422 
python esass.py 423 
python esass.py 424 
python esass.py 425 
python esass.py 426 
python esass.py 427 
python esass.py 428 
python esass.py 429 
python esass.py 430 
python esass.py 431 
python esass.py 432 
python esass.py 433 
python esass.py 434 
python esass.py 435 
python esass.py 436 
python esass.py 437 
python esass.py 438 
python esass.py 439 
python esass.py 440 
python esass.py 441 
python esass.py 442 
python esass.py 443 
python esass.py 444 
python esass.py 445 
python esass.py 446 
python esass.py 447 
python esass.py 448 
python esass.py 449 
python esass.py 450 
python esass.py 451 
python esass.py 452 
python esass.py 453 
python esass.py 454 
python esass.py 455 
python esass.py 456 
python esass.py 457 
python esass.py 458 
python esass.py 459 
python esass.py 460 
python esass.py 461 
python esass.py 462 
python esass.py 463 
python esass.py 464 
python esass.py 465 
python esass.py 466 
python esass.py 467 
python esass.py 468 
python esass.py 469 
python esass.py 470 
python esass.py 471 
python esass.py 472 
python esass.py 473 
python esass.py 474 
python esass.py 475 
python esass.py 476 
python esass.py 477 
python esass.py 478 
python esass.py 479 
python esass.py 480 
python esass.py 481 
python esass.py 482 
python esass.py 483 
python esass.py 484 
python esass.py 485 
python esass.py 486 
python esass.py 487 
python esass.py 488 
python esass.py 489 
python esass.py 490 
python esass.py 491 
python esass.py 492 
python esass.py 493 
python esass.py 494 
python esass.py 495 
python esass.py 496 
python esass.py 497 
python esass.py 498 
python esass.py 499 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 500 
python esass.py 501 
python esass.py 502 
python esass.py 503 
python esass.py 504 
python esass.py 505 
python esass.py 506 
python esass.py 507 
python esass.py 508 
python esass.py 509 
python esass.py 510 
python esass.py 511 
python esass.py 512 
python esass.py 513 
python esass.py 514 
python esass.py 515 
python esass.py 516 
python esass.py 517 
python esass.py 518 
python esass.py 519 
python esass.py 520 
python esass.py 521 
python esass.py 522 
python esass.py 523 
python esass.py 524 
python esass.py 525 
python esass.py 526 
python esass.py 527 
python esass.py 528 
python esass.py 529 
python esass.py 530 
python esass.py 531 
python esass.py 532 
python esass.py 533 
python esass.py 534 
python esass.py 535 
python esass.py 536 
python esass.py 537 
python esass.py 538 
python esass.py 539 
python esass.py 540 
python esass.py 541 
python esass.py 542 
python esass.py 543 
python esass.py 544 
python esass.py 545 
python esass.py 546 
python esass.py 547 
python esass.py 548 
python esass.py 549 
python esass.py 550 
python esass.py 551 
python esass.py 552 
python esass.py 553 
python esass.py 554 
python esass.py 555 
python esass.py 556 
python esass.py 557 
python esass.py 558 
python esass.py 559 
python esass.py 560 
python esass.py 561 
python esass.py 562 
python esass.py 563 
python esass.py 564 
python esass.py 565 
python esass.py 566 
python esass.py 567 
python esass.py 568 
python esass.py 569 
python esass.py 570 
python esass.py 571 
python esass.py 572 
python esass.py 573 
python esass.py 574 
python esass.py 575 
python esass.py 576 
python esass.py 577 
python esass.py 578 
python esass.py 579 
python esass.py 580 
python esass.py 581 
python esass.py 582 
python esass.py 583 
python esass.py 584 
python esass.py 585 
python esass.py 586 
python esass.py 587 
python esass.py 588 
python esass.py 589 
python esass.py 590 
python esass.py 591 
python esass.py 592 
python esass.py 593 
python esass.py 594 
python esass.py 595 
python esass.py 596 
python esass.py 597 
python esass.py 598 
python esass.py 599 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 600 
python esass.py 601 
python esass.py 602 
python esass.py 603 
python esass.py 604 
python esass.py 605 
python esass.py 606 
python esass.py 607 
python esass.py 608 
python esass.py 609 
python esass.py 610 
python esass.py 611 
python esass.py 612 
python esass.py 613 
python esass.py 614 
python esass.py 615 
python esass.py 616 
python esass.py 617 
python esass.py 618 
python esass.py 619 
python esass.py 620 
python esass.py 621 
python esass.py 622 
python esass.py 623 
python esass.py 624 
python esass.py 625 
python esass.py 626 
python esass.py 627 
python esass.py 628 
python esass.py 629 
python esass.py 630 
python esass.py 631 
python esass.py 632 
python esass.py 633 
python esass.py 634 
python esass.py 635 
python esass.py 636 
python esass.py 637 
python esass.py 638 
python esass.py 639 
python esass.py 640 
python esass.py 641 
python esass.py 642 
python esass.py 643 
python esass.py 644 
python esass.py 645 
python esass.py 646 
python esass.py 647 
python esass.py 648 
python esass.py 649 
python esass.py 650 
python esass.py 651 
python esass.py 652 
python esass.py 653 
python esass.py 654 
python esass.py 655 
python esass.py 656 
python esass.py 657 
python esass.py 658 
python esass.py 659 
python esass.py 660 
python esass.py 661 
python esass.py 662 
python esass.py 663 
python esass.py 664 
python esass.py 665 
python esass.py 666 
python esass.py 667 
python esass.py 668 
python esass.py 669 
python esass.py 670 
python esass.py 671 
python esass.py 672 
python esass.py 673 
python esass.py 674 
python esass.py 675 
python esass.py 676 
python esass.py 677 
python esass.py 678 
python esass.py 679 
python esass.py 680 
python esass.py 681 
python esass.py 682 
python esass.py 683 
python esass.py 684 
python esass.py 685 
python esass.py 686 
python esass.py 687 
python esass.py 688 
python esass.py 689 
python esass.py 690 
python esass.py 691 
python esass.py 692 
python esass.py 693 
python esass.py 694 
python esass.py 695 
python esass.py 696 
python esass.py 697 
python esass.py 698 
python esass.py 699 
!/bin/bash 
pyCONDA
cd $GIT_AGN_MOCK/python/sixte
source /home/erosita/sw/sass-setup.sh eSASSusers_190925 
python esass.py 700 
python esass.py 701 
python esass.py 702 
python esass.py 703 
python esass.py 704 
python esass.py 705 
python esass.py 706 
python esass.py 707 
python esass.py 708 
python esass.py 709 
python esass.py 710 
python esass.py 711 
python esass.py 712 
python esass.py 713 
python esass.py 714 
python esass.py 715 
python esass.py 716 
python esass.py 717 
python esass.py 718 
python esass.py 719 
python esass.py 720 
python esass.py 721 
python esass.py 722 
python esass.py 723 
python esass.py 724 
python esass.py 725 
python esass.py 726 
python esass.py 727 
python esass.py 728 
python esass.py 729 
python esass.py 730 
python esass.py 731 
python esass.py 732 
python esass.py 733 
python esass.py 734 
python esass.py 735 
python esass.py 736 
python esass.py 737 
python esass.py 738 
python esass.py 739 
python esass.py 740 
python esass.py 741 
python esass.py 742 
python esass.py 743 
python esass.py 744 
python esass.py 745 
python esass.py 746 
python esass.py 747 
python esass.py 748 
python esass.py 749 
python esass.py 750 
python esass.py 751 
python esass.py 752 
python esass.py 753 
python esass.py 754 
python esass.py 755 
python esass.py 756 
python esass.py 757 
python esass.py 758 
python esass.py 759 
python esass.py 760 
python esass.py 761 
python esass.py 762 
python esass.py 763 
python esass.py 764 
python esass.py 765 
python esass.py 766 
python esass.py 767 
