#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 274 
python ../pre-process-esass.py 274 
python ../esass.py 274 
