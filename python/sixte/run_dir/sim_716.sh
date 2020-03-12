#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 716 
python ../pre-process-esass.py 716 
python ../esass.py 716 
