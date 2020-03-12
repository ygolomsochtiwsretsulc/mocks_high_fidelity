#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 454 
python ../pre-process-esass.py 454 
python ../esass.py 454 
