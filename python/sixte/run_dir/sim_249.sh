#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 249 
python ../pre-process-esass.py 249 
python ../esass.py 249 
