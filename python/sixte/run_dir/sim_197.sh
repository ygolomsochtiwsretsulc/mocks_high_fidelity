#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 197 
python ../pre-process-esass.py 197 
python ../esass.py 197 
