#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 628 
python ../pre-process-esass.py 628 
python ../esass.py 628 
