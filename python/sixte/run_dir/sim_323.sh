#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 323 
python ../pre-process-esass.py 323 
python ../esass.py 323 
