#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 311 
python ../pre-process-esass.py 311 
python ../esass.py 311 
