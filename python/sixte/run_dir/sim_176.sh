#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 176 
python ../pre-process-esass.py 176 
python ../esass.py 176 
