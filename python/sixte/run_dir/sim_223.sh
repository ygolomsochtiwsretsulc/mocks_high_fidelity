#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 223 
python ../pre-process-esass.py 223 
python ../esass.py 223 
