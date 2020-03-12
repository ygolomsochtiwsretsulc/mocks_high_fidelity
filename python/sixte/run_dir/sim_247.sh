#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 247 
python ../pre-process-esass.py 247 
python ../esass.py 247 
