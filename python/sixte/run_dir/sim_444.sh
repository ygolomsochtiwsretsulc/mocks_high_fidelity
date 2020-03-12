#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 444 
python ../pre-process-esass.py 444 
python ../esass.py 444 
