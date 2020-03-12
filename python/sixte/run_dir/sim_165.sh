#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 165 
python ../pre-process-esass.py 165 
python ../esass.py 165 
