#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 381 
python ../pre-process-esass.py 381 
python ../esass.py 381 
