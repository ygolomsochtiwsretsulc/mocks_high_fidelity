#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 401 
python ../pre-process-esass.py 401 
python ../esass.py 401 
