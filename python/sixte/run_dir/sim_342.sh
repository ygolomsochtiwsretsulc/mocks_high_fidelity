#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 342 
python ../pre-process-esass.py 342 
python ../esass.py 342 
