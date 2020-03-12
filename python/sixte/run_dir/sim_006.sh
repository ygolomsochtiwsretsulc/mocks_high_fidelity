#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 006 
python ../pre-process-esass.py 006 
python ../esass.py 006 
