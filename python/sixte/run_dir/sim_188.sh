#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 188 
python ../pre-process-esass.py 188 
python ../esass.py 188 
