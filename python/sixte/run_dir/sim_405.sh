#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 405 
python ../pre-process-esass.py 405 
python ../esass.py 405 
