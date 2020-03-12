#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 756 
python ../pre-process-esass.py 756 
python ../esass.py 756 
