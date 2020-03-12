#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 135 
python ../pre-process-esass.py 135 
python ../esass.py 135 
