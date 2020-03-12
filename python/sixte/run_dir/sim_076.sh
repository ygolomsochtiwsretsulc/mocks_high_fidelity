#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 076 
python ../pre-process-esass.py 076 
python ../esass.py 076 
