#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 336 
python ../pre-process-esass.py 336 
python ../esass.py 336 
