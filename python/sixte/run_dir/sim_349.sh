#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 349 
python ../pre-process-esass.py 349 
python ../esass.py 349 
