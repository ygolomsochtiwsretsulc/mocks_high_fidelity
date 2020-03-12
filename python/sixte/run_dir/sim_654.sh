#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 654 
python ../pre-process-esass.py 654 
python ../esass.py 654 
