#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 379 
python ../pre-process-esass.py 379 
python ../esass.py 379 
