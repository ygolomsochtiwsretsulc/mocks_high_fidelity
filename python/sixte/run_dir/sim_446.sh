#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 446 
python ../pre-process-esass.py 446 
python ../esass.py 446 
