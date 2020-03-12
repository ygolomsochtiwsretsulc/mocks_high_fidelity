#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 359 
python ../pre-process-esass.py 359 
python ../esass.py 359 
