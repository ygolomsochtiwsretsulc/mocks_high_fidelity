#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 417 
python ../pre-process-esass.py 417 
python ../esass.py 417 
