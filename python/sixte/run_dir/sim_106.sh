#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 106 
python ../pre-process-esass.py 106 
python ../esass.py 106 
