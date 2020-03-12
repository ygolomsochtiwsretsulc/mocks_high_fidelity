#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 170 
python ../pre-process-esass.py 170 
python ../esass.py 170 
