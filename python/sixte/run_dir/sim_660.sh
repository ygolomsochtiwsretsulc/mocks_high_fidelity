#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 660 
python ../pre-process-esass.py 660 
python ../esass.py 660 
