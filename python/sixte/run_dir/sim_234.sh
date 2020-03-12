#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 234 
python ../pre-process-esass.py 234 
python ../esass.py 234 
