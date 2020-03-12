#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 408 
python ../pre-process-esass.py 408 
python ../esass.py 408 
