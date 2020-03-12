#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 134 
python ../pre-process-esass.py 134 
python ../esass.py 134 
