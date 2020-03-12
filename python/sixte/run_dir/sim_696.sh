#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 696 
python ../pre-process-esass.py 696 
python ../esass.py 696 
