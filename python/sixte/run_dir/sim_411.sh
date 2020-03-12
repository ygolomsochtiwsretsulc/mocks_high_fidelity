#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 411 
python ../pre-process-esass.py 411 
python ../esass.py 411 
