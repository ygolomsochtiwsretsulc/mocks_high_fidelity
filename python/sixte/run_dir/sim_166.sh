#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 166 
python ../pre-process-esass.py 166 
python ../esass.py 166 
