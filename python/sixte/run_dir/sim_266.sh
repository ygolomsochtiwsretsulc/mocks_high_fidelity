#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 266 
python ../pre-process-esass.py 266 
python ../esass.py 266 
