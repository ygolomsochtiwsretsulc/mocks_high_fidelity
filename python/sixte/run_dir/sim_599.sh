#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 599 
python ../pre-process-esass.py 599 
python ../esass.py 599 
