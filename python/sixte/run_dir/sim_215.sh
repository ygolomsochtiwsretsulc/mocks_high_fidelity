#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 215 
python ../pre-process-esass.py 215 
python ../esass.py 215 
