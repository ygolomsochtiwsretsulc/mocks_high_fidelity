#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 151 
python ../pre-process-esass.py 151 
python ../esass.py 151 
