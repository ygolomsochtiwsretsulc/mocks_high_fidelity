#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 498 
python ../pre-process-esass.py 498 
python ../esass.py 498 
