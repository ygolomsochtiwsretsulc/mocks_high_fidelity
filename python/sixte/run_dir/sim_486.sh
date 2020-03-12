#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 486 
python ../pre-process-esass.py 486 
python ../esass.py 486 
