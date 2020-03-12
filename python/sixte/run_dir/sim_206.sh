#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 206 
python ../pre-process-esass.py 206 
python ../esass.py 206 
