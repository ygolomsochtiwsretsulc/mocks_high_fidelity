#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 300 
python ../pre-process-esass.py 300 
python ../esass.py 300 
