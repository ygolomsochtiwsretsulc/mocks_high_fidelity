#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 580 
python ../pre-process-esass.py 580 
python ../esass.py 580 
