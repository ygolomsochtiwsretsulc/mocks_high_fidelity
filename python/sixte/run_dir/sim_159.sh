#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 159 
python ../pre-process-esass.py 159 
python ../esass.py 159 
