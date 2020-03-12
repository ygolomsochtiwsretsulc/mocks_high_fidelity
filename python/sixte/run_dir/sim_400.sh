#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 400 
python ../pre-process-esass.py 400 
python ../esass.py 400 
