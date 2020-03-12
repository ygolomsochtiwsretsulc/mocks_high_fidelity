#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 108 
python ../pre-process-esass.py 108 
python ../esass.py 108 
