#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 097 
python ../pre-process-esass.py 097 
python ../esass.py 097 
