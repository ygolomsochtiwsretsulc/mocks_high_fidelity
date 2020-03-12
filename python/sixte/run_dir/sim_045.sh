#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 045 
python ../pre-process-esass.py 045 
python ../esass.py 045 
