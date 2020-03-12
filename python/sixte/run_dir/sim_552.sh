#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 552 
python ../pre-process-esass.py 552 
python ../esass.py 552 
