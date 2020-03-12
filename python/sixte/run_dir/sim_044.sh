#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 044 
python ../pre-process-esass.py 044 
python ../esass.py 044 
