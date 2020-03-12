#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 546 
python ../pre-process-esass.py 546 
python ../esass.py 546 
