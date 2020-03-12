#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 240 
python ../pre-process-esass.py 240 
python ../esass.py 240 
