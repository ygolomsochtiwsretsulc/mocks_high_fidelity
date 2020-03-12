#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 265 
python ../pre-process-esass.py 265 
python ../esass.py 265 
