#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 519 
python ../pre-process-esass.py 519 
python ../esass.py 519 
