#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 081 
python ../pre-process-esass.py 081 
python ../esass.py 081 
