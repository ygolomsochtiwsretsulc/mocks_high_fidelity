#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 456 
python ../pre-process-esass.py 456 
python ../esass.py 456 
