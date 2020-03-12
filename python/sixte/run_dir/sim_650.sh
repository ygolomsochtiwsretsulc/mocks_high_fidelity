#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 650 
python ../pre-process-esass.py 650 
python ../esass.py 650 
