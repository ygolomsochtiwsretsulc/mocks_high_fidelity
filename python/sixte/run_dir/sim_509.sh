#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 509 
python ../pre-process-esass.py 509 
python ../esass.py 509 
