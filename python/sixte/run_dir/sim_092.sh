#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 092 
python ../pre-process-esass.py 092 
python ../esass.py 092 
