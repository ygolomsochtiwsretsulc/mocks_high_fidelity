#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 632 
python ../pre-process-esass.py 632 
python ../esass.py 632 
