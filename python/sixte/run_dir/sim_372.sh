#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 372 
python ../pre-process-esass.py 372 
python ../esass.py 372 
