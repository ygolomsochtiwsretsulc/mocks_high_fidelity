#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 250 
python ../pre-process-esass.py 250 
python ../esass.py 250 
