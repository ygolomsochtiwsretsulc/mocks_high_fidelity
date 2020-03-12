#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 028 
python ../pre-process-esass.py 028 
python ../esass.py 028 
