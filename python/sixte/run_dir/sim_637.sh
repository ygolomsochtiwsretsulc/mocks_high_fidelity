#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 637 
python ../pre-process-esass.py 637 
python ../esass.py 637 
