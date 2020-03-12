#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 100 
python ../pre-process-esass.py 100 
python ../esass.py 100 
