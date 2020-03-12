#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 730 
python ../pre-process-esass.py 730 
python ../esass.py 730 
