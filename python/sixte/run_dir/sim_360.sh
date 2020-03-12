#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 360 
python ../pre-process-esass.py 360 
python ../esass.py 360 
