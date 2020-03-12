#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 207 
python ../pre-process-esass.py 207 
python ../esass.py 207 
