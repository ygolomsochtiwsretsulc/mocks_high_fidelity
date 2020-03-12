#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 540 
python ../pre-process-esass.py 540 
python ../esass.py 540 
