#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 707 
python ../pre-process-esass.py 707 
python ../esass.py 707 
