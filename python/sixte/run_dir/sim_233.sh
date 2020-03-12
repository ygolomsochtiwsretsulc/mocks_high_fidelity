#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 233 
python ../pre-process-esass.py 233 
python ../esass.py 233 
