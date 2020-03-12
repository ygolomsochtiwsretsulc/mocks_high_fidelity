#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 180 
python ../pre-process-esass.py 180 
python ../esass.py 180 
