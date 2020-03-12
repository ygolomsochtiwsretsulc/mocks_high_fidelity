#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 370 
python ../pre-process-esass.py 370 
python ../esass.py 370 
