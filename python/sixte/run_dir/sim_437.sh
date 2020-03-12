#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 437 
python ../pre-process-esass.py 437 
python ../esass.py 437 
