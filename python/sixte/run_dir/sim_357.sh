#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 357 
python ../pre-process-esass.py 357 
python ../esass.py 357 
