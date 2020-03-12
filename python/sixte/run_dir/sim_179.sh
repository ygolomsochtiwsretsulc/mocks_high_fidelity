#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 179 
python ../pre-process-esass.py 179 
python ../esass.py 179 
