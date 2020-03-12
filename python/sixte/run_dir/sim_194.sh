#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 194 
python ../pre-process-esass.py 194 
python ../esass.py 194 
