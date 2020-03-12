#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 365 
python ../pre-process-esass.py 365 
python ../esass.py 365 
