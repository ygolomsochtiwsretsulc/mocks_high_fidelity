#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 228 
python ../pre-process-esass.py 228 
python ../esass.py 228 
