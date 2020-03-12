#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 474 
python ../pre-process-esass.py 474 
python ../esass.py 474 
