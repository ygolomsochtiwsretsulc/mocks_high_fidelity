#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 421 
python ../pre-process-esass.py 421 
python ../esass.py 421 
