#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 558 
python ../pre-process-esass.py 558 
python ../esass.py 558 
