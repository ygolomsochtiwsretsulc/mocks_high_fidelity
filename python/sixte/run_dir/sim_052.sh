#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 052 
python ../pre-process-esass.py 052 
python ../esass.py 052 
