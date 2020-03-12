#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 190 
python ../pre-process-esass.py 190 
python ../esass.py 190 
