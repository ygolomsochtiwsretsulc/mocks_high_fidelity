#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 137 
python ../pre-process-esass.py 137 
python ../esass.py 137 
