#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 039 
python ../pre-process-esass.py 039 
python ../esass.py 039 
