#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 701 
python ../pre-process-esass.py 701 
python ../esass.py 701 
