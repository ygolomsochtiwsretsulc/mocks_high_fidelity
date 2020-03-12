#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 220 
python ../pre-process-esass.py 220 
python ../esass.py 220 
