#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 369 
python ../pre-process-esass.py 369 
python ../esass.py 369 
