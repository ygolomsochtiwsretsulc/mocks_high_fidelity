#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 625 
python ../pre-process-esass.py 625 
python ../esass.py 625 
