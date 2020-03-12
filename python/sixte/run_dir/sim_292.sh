#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 292 
python ../pre-process-esass.py 292 
python ../esass.py 292 
