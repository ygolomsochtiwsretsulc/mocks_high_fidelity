#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 221 
python ../pre-process-esass.py 221 
python ../esass.py 221 
