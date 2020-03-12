#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 380 
python ../pre-process-esass.py 380 
python ../esass.py 380 
