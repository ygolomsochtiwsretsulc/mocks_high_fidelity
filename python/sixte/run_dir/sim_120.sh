#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 120 
python ../pre-process-esass.py 120 
python ../esass.py 120 
