#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 697 
python ../pre-process-esass.py 697 
python ../esass.py 697 
