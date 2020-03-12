#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 528 
python ../pre-process-esass.py 528 
python ../esass.py 528 
