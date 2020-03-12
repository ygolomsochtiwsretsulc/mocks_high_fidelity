#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 301 
python ../pre-process-esass.py 301 
python ../esass.py 301 
