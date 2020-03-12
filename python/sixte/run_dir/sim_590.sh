#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 590 
python ../pre-process-esass.py 590 
python ../esass.py 590 
