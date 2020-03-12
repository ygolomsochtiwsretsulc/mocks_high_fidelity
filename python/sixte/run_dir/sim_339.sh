#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 339 
python ../pre-process-esass.py 339 
python ../esass.py 339 
