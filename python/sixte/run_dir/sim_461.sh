#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 461 
python ../pre-process-esass.py 461 
python ../esass.py 461 
