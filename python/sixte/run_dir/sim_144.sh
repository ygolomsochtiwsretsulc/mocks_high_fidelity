#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 144 
python ../pre-process-esass.py 144 
python ../esass.py 144 
