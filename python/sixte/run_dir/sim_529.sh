#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 529 
python ../pre-process-esass.py 529 
python ../esass.py 529 
