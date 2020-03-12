#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 341 
python ../pre-process-esass.py 341 
python ../esass.py 341 
