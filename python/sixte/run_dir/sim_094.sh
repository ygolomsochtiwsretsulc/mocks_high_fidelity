#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 094 
python ../pre-process-esass.py 094 
python ../esass.py 094 
