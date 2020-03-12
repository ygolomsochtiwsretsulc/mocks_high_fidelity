#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 531 
python ../pre-process-esass.py 531 
python ../esass.py 531 
