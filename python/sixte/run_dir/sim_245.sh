#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 245 
python ../pre-process-esass.py 245 
python ../esass.py 245 
