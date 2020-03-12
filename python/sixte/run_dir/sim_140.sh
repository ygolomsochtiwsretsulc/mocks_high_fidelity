#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 140 
python ../pre-process-esass.py 140 
python ../esass.py 140 
