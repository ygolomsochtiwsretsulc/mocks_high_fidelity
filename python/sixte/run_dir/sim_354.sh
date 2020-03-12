#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 354 
python ../pre-process-esass.py 354 
python ../esass.py 354 
