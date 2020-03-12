#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 218 
python ../pre-process-esass.py 218 
python ../esass.py 218 
