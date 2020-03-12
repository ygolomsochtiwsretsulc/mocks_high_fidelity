#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 231 
python ../pre-process-esass.py 231 
python ../esass.py 231 
