#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 212 
python ../pre-process-esass.py 212 
python ../esass.py 212 
