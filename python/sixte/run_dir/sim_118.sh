#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 118 
python ../pre-process-esass.py 118 
python ../esass.py 118 
