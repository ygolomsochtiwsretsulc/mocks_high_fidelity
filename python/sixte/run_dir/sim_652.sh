#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 652 
python ../pre-process-esass.py 652 
python ../esass.py 652 
