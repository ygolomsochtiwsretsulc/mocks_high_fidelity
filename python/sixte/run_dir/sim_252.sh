#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 252 
python ../pre-process-esass.py 252 
python ../esass.py 252 
