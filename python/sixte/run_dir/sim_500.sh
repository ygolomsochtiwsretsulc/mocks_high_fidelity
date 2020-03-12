#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 500 
python ../pre-process-esass.py 500 
python ../esass.py 500 
