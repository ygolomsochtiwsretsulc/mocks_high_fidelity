#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 239 
python ../pre-process-esass.py 239 
python ../esass.py 239 
