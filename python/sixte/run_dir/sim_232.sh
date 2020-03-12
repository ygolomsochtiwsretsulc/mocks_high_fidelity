#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 232 
python ../pre-process-esass.py 232 
python ../esass.py 232 
