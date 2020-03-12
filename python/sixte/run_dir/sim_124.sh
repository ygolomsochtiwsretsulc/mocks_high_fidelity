#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 124 
python ../pre-process-esass.py 124 
python ../esass.py 124 
