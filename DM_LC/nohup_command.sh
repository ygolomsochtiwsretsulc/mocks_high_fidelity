#!/bin/bash

# RUN ON DS43
pyCONDA 

python 000_geometry.py 1000.0 MD10
python 000_geometry.py 400.0  MD04
python 000_geometry.py 4000.0 MD40	
python 000_geometry.py 1000.0 UNIT_fA1_DIR 
python 000_geometry.py 1000.0 UNIT_fA2_DIR 
python 000_geometry.py 1000.0 UNIT_fA1i_DIR
python 000_geometry.py 1000.0 UNIT_fA2i_DIR

python 001_process_hlists_WRITE_commands.py
# copy paste outputs in the md??_r000_all.sh
# Runs started on Dec. 12
# MD40, MDR1, MD25 and UNIT_fA2i_DIR still downloading at MPCDF, soo will need to be rerun
# ALL z=0 UNIT snapshots re - downloading on ds54 
# need to rerun them once completed
nohup sh md04_r000_all.sh > md04_r000_all.log & # onoging ds43
nohup sh md10_r000_all.sh > md10_r000_all.log & # ongoing ds43
nohup sh md40_r000_all.sh > md40_r000_all.log & # ongoing ds52
nohup sh md40_r000_all_2.sh > md40_r000_all_2.log & # ongoing ds52
nohup sh unit_fA1i_r000_all.sh  > fA1i_r000_all.log & # ds43 onoging
nohup sh unit_fA1_r000_all.sh   > fA1_r000_all.log & # ds43 onoging
nohup sh unit_fA2_r000_all.sh   > fA2_r000_all.log & # ds43 onoging

