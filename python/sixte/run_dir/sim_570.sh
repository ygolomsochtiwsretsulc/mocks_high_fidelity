#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 570 
python ../pre-process-esass.py 570 
python ../esass.py 570 
