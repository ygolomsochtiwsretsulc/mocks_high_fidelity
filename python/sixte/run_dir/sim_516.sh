#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 516 
python ../pre-process-esass.py 516 
python ../esass.py 516 
