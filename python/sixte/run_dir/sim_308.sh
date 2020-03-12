#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 308 
python ../pre-process-esass.py 308 
python ../esass.py 308 
