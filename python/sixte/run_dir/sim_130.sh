#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 130 
python ../pre-process-esass.py 130 
python ../esass.py 130 
