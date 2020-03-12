#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 574 
python ../pre-process-esass.py 574 
python ../esass.py 574 
