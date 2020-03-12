#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 695 
python ../pre-process-esass.py 695 
python ../esass.py 695 
