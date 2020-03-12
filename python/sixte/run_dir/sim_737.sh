#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 737 
python ../pre-process-esass.py 737 
python ../esass.py 737 
