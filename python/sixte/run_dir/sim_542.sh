#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 542 
python ../pre-process-esass.py 542 
python ../esass.py 542 
