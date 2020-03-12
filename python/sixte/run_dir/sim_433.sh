#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 433 
python ../pre-process-esass.py 433 
python ../esass.py 433 
