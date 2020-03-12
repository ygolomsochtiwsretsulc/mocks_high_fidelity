#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 316 
python ../pre-process-esass.py 316 
python ../esass.py 316 
