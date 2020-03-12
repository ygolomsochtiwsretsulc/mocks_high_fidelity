#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 157 
python ../pre-process-esass.py 157 
python ../esass.py 157 
