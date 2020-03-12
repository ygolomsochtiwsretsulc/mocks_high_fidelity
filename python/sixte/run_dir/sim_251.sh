#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 251 
python ../pre-process-esass.py 251 
python ../esass.py 251 
