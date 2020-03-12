#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 743 
python ../pre-process-esass.py 743 
python ../esass.py 743 
