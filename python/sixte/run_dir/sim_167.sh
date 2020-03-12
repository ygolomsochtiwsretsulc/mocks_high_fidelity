#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 167 
python ../pre-process-esass.py 167 
python ../esass.py 167 
