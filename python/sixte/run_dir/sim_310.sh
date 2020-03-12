#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 310 
python ../pre-process-esass.py 310 
python ../esass.py 310 
