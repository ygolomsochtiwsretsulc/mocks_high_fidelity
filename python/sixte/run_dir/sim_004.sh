#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 004 
python ../pre-process-esass.py 004 
python ../esass.py 004 
