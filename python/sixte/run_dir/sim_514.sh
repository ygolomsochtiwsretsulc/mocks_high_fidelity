#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 514 
python ../pre-process-esass.py 514 
python ../esass.py 514 
