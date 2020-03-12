#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 018 
python ../pre-process-esass.py 018 
python ../esass.py 018 
