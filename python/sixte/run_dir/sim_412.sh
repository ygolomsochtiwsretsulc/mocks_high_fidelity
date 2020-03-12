#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 412 
python ../pre-process-esass.py 412 
python ../esass.py 412 
