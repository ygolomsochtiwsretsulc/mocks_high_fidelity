#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 399 
python ../pre-process-esass.py 399 
python ../esass.py 399 
