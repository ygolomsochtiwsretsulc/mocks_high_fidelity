#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 351 
python ../pre-process-esass.py 351 
python ../esass.py 351 
