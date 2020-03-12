#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 295 
python ../pre-process-esass.py 295 
python ../esass.py 295 
