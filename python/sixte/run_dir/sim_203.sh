#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 203 
python ../pre-process-esass.py 203 
python ../esass.py 203 
