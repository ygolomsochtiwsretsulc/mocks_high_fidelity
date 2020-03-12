#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 199 
python ../pre-process-esass.py 199 
python ../esass.py 199 
