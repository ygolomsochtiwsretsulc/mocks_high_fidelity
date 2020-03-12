#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 297 
python ../pre-process-esass.py 297 
python ../esass.py 297 
