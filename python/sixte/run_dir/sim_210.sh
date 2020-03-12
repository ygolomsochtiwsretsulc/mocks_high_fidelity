#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 210 
python ../pre-process-esass.py 210 
python ../esass.py 210 
