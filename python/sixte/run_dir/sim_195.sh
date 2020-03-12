#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 195 
python ../pre-process-esass.py 195 
python ../esass.py 195 
