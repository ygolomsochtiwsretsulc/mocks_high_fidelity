#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 314 
python ../pre-process-esass.py 314 
python ../esass.py 314 
