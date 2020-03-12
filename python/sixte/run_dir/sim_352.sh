#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 352 
python ../pre-process-esass.py 352 
python ../esass.py 352 
