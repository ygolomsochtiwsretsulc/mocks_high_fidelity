#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 608 
python ../pre-process-esass.py 608 
python ../esass.py 608 
