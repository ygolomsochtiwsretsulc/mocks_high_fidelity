#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 288 
python ../pre-process-esass.py 288 
python ../esass.py 288 
