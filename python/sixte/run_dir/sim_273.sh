#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 273 
python ../pre-process-esass.py 273 
python ../esass.py 273 
