#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 277 
python ../pre-process-esass.py 277 
python ../esass.py 277 
