#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 294 
python ../pre-process-esass.py 294 
python ../esass.py 294 
