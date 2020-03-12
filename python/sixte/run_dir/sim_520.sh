#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 520 
python ../pre-process-esass.py 520 
python ../esass.py 520 
