#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 101 
python ../pre-process-esass.py 101 
python ../esass.py 101 
