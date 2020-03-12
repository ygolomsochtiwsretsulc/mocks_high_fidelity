#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 325 
python ../pre-process-esass.py 325 
python ../esass.py 325 
