#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 112 
python ../pre-process-esass.py 112 
python ../esass.py 112 
