#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 410 
python ../pre-process-esass.py 410 
python ../esass.py 410 
