#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 332 
python ../pre-process-esass.py 332 
python ../esass.py 332 
