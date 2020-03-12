#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 263 
python ../pre-process-esass.py 263 
python ../esass.py 263 
