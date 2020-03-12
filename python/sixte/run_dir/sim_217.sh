#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 217 
python ../pre-process-esass.py 217 
python ../esass.py 217 
