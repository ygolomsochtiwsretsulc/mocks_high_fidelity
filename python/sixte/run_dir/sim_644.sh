#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 644 
python ../pre-process-esass.py 644 
python ../esass.py 644 
