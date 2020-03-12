#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 299 
python ../pre-process-esass.py 299 
python ../esass.py 299 
