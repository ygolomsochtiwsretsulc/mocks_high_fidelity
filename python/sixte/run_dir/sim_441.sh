#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 441 
python ../pre-process-esass.py 441 
python ../esass.py 441 
