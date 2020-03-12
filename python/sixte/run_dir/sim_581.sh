#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 581 
python ../pre-process-esass.py 581 
python ../esass.py 581 
