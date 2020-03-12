#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 279 
python ../pre-process-esass.py 279 
python ../esass.py 279 
