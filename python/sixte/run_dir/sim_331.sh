#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 331 
python ../pre-process-esass.py 331 
python ../esass.py 331 
