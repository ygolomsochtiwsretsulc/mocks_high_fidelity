#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 504 
python ../pre-process-esass.py 504 
python ../esass.py 504 
