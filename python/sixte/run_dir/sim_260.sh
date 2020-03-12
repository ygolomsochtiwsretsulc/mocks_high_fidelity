#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 260 
python ../pre-process-esass.py 260 
python ../esass.py 260 
