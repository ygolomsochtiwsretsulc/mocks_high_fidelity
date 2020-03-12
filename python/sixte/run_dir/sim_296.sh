#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 296 
python ../pre-process-esass.py 296 
python ../esass.py 296 
