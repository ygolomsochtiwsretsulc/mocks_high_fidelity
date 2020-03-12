#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 448 
python ../pre-process-esass.py 448 
python ../esass.py 448 
