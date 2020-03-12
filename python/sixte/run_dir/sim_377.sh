#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 377 
python ../pre-process-esass.py 377 
python ../esass.py 377 
