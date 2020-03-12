#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 521 
python ../pre-process-esass.py 521 
python ../esass.py 521 
