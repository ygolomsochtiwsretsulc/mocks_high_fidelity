#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 205 
python ../pre-process-esass.py 205 
python ../esass.py 205 
