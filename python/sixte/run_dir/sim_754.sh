#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 754 
python ../pre-process-esass.py 754 
python ../esass.py 754 
