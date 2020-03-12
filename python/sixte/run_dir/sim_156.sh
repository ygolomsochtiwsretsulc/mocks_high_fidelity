#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 156 
python ../pre-process-esass.py 156 
python ../esass.py 156 
