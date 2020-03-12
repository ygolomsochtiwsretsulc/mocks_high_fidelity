#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 703 
python ../pre-process-esass.py 703 
python ../esass.py 703 
