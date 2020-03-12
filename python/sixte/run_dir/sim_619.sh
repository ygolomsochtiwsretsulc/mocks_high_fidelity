#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 619 
python ../pre-process-esass.py 619 
python ../esass.py 619 
