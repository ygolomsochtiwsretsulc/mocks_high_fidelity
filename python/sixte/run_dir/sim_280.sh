#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 280 
python ../pre-process-esass.py 280 
python ../esass.py 280 
