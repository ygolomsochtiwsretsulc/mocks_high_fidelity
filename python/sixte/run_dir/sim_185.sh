#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 185 
python ../pre-process-esass.py 185 
python ../esass.py 185 
