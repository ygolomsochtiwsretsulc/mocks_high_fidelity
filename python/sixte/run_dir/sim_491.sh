#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 491 
python ../pre-process-esass.py 491 
python ../esass.py 491 
