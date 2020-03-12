#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 524 
python ../pre-process-esass.py 524 
python ../esass.py 524 
