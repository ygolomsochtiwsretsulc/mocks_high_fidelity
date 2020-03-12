#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 196 
python ../pre-process-esass.py 196 
python ../esass.py 196 
