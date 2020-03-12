#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 258 
python ../pre-process-esass.py 258 
python ../esass.py 258 
