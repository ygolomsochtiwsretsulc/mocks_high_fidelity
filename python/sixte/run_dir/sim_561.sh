#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 561 
python ../pre-process-esass.py 561 
python ../esass.py 561 
