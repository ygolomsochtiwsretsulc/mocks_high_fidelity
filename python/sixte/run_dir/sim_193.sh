#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 193 
python ../pre-process-esass.py 193 
python ../esass.py 193 
