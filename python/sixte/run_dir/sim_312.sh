#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 312 
python ../pre-process-esass.py 312 
python ../esass.py 312 
