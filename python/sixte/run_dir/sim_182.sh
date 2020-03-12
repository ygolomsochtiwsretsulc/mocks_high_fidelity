#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 182 
python ../pre-process-esass.py 182 
python ../esass.py 182 
