#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 329 
python ../pre-process-esass.py 329 
python ../esass.py 329 
