#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 677 
python ../pre-process-esass.py 677 
python ../esass.py 677 
