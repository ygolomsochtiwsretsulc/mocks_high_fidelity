#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 714 
python ../pre-process-esass.py 714 
python ../esass.py 714 
