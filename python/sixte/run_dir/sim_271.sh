#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 271 
python ../pre-process-esass.py 271 
python ../esass.py 271 
