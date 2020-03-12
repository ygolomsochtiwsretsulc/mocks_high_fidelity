#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 404 
python ../pre-process-esass.py 404 
python ../esass.py 404 
