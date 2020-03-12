#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 661 
python ../pre-process-esass.py 661 
python ../esass.py 661 
