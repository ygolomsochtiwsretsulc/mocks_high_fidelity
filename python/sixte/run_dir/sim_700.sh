#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 700 
python ../pre-process-esass.py 700 
python ../esass.py 700 
