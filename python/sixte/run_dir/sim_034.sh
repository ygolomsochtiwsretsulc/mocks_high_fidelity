#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 034 
python ../pre-process-esass.py 034 
python ../esass.py 034 
