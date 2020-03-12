#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 059 
python ../pre-process-esass.py 059 
python ../esass.py 059 
