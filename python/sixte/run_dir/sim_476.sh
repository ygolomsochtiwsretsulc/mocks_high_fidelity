#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 476 
python ../pre-process-esass.py 476 
python ../esass.py 476 
