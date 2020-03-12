#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 458 
python ../pre-process-esass.py 458 
python ../esass.py 458 
