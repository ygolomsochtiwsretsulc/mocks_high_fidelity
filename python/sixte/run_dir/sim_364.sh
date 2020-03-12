#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 364 
python ../pre-process-esass.py 364 
python ../esass.py 364 
