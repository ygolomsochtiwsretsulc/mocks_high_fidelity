#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 161 
python ../pre-process-esass.py 161 
python ../esass.py 161 
