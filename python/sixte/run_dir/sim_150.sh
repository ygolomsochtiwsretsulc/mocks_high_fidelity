#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 150 
python ../pre-process-esass.py 150 
python ../esass.py 150 
