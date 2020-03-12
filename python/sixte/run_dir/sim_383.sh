#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 383 
python ../pre-process-esass.py 383 
python ../esass.py 383 
