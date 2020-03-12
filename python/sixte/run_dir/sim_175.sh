#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 175 
python ../pre-process-esass.py 175 
python ../esass.py 175 
