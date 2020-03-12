#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 327 
python ../pre-process-esass.py 327 
python ../esass.py 327 
