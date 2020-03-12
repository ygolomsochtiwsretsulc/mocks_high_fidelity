#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 202 
python ../pre-process-esass.py 202 
python ../esass.py 202 
