#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 069 
python ../pre-process-esass.py 069 
python ../esass.py 069 
