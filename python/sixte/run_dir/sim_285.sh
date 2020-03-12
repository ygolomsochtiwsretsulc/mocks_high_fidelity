#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 285 
python ../pre-process-esass.py 285 
python ../esass.py 285 
