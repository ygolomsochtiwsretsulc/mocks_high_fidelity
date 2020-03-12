#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 398 
python ../pre-process-esass.py 398 
python ../esass.py 398 
