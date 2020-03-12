#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 738 
python ../pre-process-esass.py 738 
python ../esass.py 738 
