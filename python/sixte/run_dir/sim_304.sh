#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 304 
python ../pre-process-esass.py 304 
python ../esass.py 304 
