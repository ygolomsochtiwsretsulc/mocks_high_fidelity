#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 162 
python ../pre-process-esass.py 162 
python ../esass.py 162 
