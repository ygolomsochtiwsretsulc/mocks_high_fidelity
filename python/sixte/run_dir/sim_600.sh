#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 600 
python ../pre-process-esass.py 600 
python ../esass.py 600 
