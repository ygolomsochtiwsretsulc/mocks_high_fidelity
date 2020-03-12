#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 407 
python ../pre-process-esass.py 407 
python ../esass.py 407 
