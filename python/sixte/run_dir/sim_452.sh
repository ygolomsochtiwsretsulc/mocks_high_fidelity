#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 452 
python ../pre-process-esass.py 452 
python ../esass.py 452 
