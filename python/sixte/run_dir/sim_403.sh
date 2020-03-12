#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 403 
python ../pre-process-esass.py 403 
python ../esass.py 403 
