#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 429 
python ../pre-process-esass.py 429 
python ../esass.py 429 
