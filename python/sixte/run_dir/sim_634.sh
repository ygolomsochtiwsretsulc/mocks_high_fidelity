#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 634 
python ../pre-process-esass.py 634 
python ../esass.py 634 
