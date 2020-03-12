#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 119 
python ../pre-process-esass.py 119 
python ../esass.py 119 
