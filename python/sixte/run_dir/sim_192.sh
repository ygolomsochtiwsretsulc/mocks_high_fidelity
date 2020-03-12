#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 192 
python ../pre-process-esass.py 192 
python ../esass.py 192 
