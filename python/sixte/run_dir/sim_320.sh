#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 320 
python ../pre-process-esass.py 320 
python ../esass.py 320 
