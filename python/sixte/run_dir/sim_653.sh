#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 653 
python ../pre-process-esass.py 653 
python ../esass.py 653 
