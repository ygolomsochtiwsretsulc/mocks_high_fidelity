#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 344 
python ../pre-process-esass.py 344 
python ../esass.py 344 
