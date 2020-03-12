#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 012 
python ../pre-process-esass.py 012 
python ../esass.py 012 
