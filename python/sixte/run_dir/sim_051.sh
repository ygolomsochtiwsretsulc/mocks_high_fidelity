#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 051 
python ../pre-process-esass.py 051 
python ../esass.py 051 
