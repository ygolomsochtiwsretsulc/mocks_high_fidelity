#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 363 
python ../pre-process-esass.py 363 
python ../esass.py 363 
