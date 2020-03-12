#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 554 
python ../pre-process-esass.py 554 
python ../esass.py 554 
