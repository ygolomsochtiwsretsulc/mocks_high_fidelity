#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 264 
python ../pre-process-esass.py 264 
python ../esass.py 264 
