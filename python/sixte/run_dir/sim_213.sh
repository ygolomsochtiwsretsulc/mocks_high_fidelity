#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 213 
python ../pre-process-esass.py 213 
python ../esass.py 213 
