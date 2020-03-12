#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 001 
python ../pre-process-esass.py 001 
python ../esass.py 001 
