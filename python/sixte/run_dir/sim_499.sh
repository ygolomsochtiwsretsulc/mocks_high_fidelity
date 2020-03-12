#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 499 
python ../pre-process-esass.py 499 
python ../esass.py 499 
