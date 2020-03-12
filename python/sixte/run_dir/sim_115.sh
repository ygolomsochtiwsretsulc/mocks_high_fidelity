#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 115 
python ../pre-process-esass.py 115 
python ../esass.py 115 
