#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 481 
python ../pre-process-esass.py 481 
python ../esass.py 481 
