#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 084 
python ../pre-process-esass.py 084 
python ../esass.py 084 
