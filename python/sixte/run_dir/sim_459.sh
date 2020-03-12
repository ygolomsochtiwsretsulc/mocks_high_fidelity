#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 459 
python ../pre-process-esass.py 459 
python ../esass.py 459 
