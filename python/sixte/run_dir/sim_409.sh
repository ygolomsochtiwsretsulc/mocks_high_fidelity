#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 409 
python ../pre-process-esass.py 409 
python ../esass.py 409 
