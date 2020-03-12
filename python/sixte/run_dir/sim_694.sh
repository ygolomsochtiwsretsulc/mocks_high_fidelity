#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 694 
python ../pre-process-esass.py 694 
python ../esass.py 694 
