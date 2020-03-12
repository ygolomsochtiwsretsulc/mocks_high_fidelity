#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 390 
python ../pre-process-esass.py 390 
python ../esass.py 390 
