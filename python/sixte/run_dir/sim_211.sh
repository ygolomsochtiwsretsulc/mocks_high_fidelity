#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 211 
python ../pre-process-esass.py 211 
python ../esass.py 211 
