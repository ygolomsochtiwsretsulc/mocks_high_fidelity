#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 436 
python ../pre-process-esass.py 436 
python ../esass.py 436 
