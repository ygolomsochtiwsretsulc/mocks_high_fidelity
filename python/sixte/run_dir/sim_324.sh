#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 324 
python ../pre-process-esass.py 324 
python ../esass.py 324 
