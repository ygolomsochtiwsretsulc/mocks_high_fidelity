#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 385 
python ../pre-process-esass.py 385 
python ../esass.py 385 
