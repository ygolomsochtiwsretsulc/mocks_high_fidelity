#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 348 
python ../pre-process-esass.py 348 
python ../esass.py 348 
