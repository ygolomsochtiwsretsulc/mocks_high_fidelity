#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 645 
python ../pre-process-esass.py 645 
python ../esass.py 645 
