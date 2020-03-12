#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 604 
python ../pre-process-esass.py 604 
python ../esass.py 604 
