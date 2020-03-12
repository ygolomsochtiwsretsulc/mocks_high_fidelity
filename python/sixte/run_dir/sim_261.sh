#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 261 
python ../pre-process-esass.py 261 
python ../esass.py 261 
