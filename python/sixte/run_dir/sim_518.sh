#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 518 
python ../pre-process-esass.py 518 
python ../esass.py 518 
