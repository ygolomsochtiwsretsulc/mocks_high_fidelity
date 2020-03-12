#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 375 
python ../pre-process-esass.py 375 
python ../esass.py 375 
