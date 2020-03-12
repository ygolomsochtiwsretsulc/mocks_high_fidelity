#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 328 
python ../pre-process-esass.py 328 
python ../esass.py 328 
