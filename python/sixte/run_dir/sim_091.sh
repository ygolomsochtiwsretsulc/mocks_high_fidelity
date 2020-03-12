#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 091 
python ../pre-process-esass.py 091 
python ../esass.py 091 
