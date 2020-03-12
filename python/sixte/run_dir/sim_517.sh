#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 517 
python ../pre-process-esass.py 517 
python ../esass.py 517 
