#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 469 
python ../pre-process-esass.py 469 
python ../esass.py 469 
