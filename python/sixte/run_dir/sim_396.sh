#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 396 
python ../pre-process-esass.py 396 
python ../esass.py 396 
