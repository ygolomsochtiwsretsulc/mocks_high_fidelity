#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 537 
python ../pre-process-esass.py 537 
python ../esass.py 537 
