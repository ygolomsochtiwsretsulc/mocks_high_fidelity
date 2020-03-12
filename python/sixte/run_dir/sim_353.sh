#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 353 
python ../pre-process-esass.py 353 
python ../esass.py 353 
