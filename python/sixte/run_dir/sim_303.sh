#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 303 
python ../pre-process-esass.py 303 
python ../esass.py 303 
