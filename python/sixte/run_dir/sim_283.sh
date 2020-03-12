#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 283 
python ../pre-process-esass.py 283 
python ../esass.py 283 
