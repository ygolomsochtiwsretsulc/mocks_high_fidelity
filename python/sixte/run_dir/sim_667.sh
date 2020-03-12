#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 667 
python ../pre-process-esass.py 667 
python ../esass.py 667 
