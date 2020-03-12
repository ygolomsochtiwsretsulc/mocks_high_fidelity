#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 525 
python ../pre-process-esass.py 525 
python ../esass.py 525 
