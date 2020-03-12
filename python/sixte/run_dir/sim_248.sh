#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 248 
python ../pre-process-esass.py 248 
python ../esass.py 248 
