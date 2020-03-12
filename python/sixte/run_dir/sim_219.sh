#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 219 
python ../pre-process-esass.py 219 
python ../esass.py 219 
