#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 145 
python ../pre-process-esass.py 145 
python ../esass.py 145 
