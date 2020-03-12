#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 259 
python ../pre-process-esass.py 259 
python ../esass.py 259 
