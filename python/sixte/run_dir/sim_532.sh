#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 532 
python ../pre-process-esass.py 532 
python ../esass.py 532 
