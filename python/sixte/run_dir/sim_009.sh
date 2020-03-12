#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 009 
python ../pre-process-esass.py 009 
python ../esass.py 009 
