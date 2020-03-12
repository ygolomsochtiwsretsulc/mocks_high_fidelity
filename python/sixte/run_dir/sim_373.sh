#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 373 
python ../pre-process-esass.py 373 
python ../esass.py 373 
