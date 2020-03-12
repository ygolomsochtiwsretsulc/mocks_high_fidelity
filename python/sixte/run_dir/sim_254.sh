#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 254 
python ../pre-process-esass.py 254 
python ../esass.py 254 
