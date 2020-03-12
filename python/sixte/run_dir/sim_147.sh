#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 147 
python ../pre-process-esass.py 147 
python ../esass.py 147 
