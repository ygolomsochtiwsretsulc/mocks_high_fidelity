#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 305 
python ../pre-process-esass.py 305 
python ../esass.py 305 
