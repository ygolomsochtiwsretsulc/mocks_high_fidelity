#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 318 
python ../pre-process-esass.py 318 
python ../esass.py 318 
