#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 343 
python ../pre-process-esass.py 343 
python ../esass.py 343 
