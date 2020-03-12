#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 128 
python ../pre-process-esass.py 128 
python ../esass.py 128 
