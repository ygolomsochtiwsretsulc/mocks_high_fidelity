#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 321 
python ../pre-process-esass.py 321 
python ../esass.py 321 
