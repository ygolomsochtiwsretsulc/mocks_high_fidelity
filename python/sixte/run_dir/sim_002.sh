#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 002 
python ../pre-process-esass.py 002 
python ../esass.py 002 
