#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 262 
python ../pre-process-esass.py 262 
python ../esass.py 262 
