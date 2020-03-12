#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 340 
python ../pre-process-esass.py 340 
python ../esass.py 340 
