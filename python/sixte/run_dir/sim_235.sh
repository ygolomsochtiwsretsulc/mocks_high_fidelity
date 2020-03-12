#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 235 
python ../pre-process-esass.py 235 
python ../esass.py 235 
