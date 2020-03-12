#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 736 
python ../pre-process-esass.py 736 
python ../esass.py 736 
