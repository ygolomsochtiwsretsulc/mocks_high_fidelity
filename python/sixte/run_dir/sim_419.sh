#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 419 
python ../pre-process-esass.py 419 
python ../esass.py 419 
