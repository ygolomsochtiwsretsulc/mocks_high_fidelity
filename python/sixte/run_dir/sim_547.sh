#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 547 
python ../pre-process-esass.py 547 
python ../esass.py 547 
