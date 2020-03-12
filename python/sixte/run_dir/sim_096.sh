#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 096 
python ../pre-process-esass.py 096 
python ../esass.py 096 
