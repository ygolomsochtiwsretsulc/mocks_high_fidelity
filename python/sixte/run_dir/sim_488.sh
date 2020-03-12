#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 488 
python ../pre-process-esass.py 488 
python ../esass.py 488 
