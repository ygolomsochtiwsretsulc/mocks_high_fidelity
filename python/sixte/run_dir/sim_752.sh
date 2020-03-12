#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 752 
python ../pre-process-esass.py 752 
python ../esass.py 752 
