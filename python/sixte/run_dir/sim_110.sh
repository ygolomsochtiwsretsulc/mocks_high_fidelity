#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 110 
python ../pre-process-esass.py 110 
python ../esass.py 110 
