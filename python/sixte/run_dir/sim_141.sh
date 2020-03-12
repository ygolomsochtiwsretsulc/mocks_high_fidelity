#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 141 
python ../pre-process-esass.py 141 
python ../esass.py 141 
