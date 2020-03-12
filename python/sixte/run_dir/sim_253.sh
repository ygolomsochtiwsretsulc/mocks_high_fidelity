#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 253 
python ../pre-process-esass.py 253 
python ../esass.py 253 
