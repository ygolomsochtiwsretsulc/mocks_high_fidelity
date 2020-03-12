#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 055 
python ../pre-process-esass.py 055 
python ../esass.py 055 
