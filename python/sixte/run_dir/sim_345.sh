#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 345 
python ../pre-process-esass.py 345 
python ../esass.py 345 
