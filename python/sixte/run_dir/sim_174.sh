#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 174 
python ../pre-process-esass.py 174 
python ../esass.py 174 
