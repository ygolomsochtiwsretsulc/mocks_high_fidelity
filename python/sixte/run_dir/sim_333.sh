#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 333 
python ../pre-process-esass.py 333 
python ../esass.py 333 
