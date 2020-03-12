#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 032 
python ../pre-process-esass.py 032 
python ../esass.py 032 
