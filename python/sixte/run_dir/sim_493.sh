#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 493 
python ../pre-process-esass.py 493 
python ../esass.py 493 
