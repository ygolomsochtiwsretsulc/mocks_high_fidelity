#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 281 
python ../pre-process-esass.py 281 
python ../esass.py 281 
