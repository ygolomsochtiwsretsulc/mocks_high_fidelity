#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 536 
python ../pre-process-esass.py 536 
python ../esass.py 536 
