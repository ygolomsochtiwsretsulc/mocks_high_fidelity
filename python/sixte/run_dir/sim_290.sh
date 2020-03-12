#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 290 
python ../pre-process-esass.py 290 
python ../esass.py 290 
