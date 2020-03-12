#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 395 
python ../pre-process-esass.py 395 
python ../esass.py 395 
