#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 389 
python ../pre-process-esass.py 389 
python ../esass.py 389 
