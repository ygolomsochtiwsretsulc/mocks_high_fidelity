#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 125 
python ../pre-process-esass.py 125 
python ../esass.py 125 
