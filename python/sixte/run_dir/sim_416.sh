#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 416 
python ../pre-process-esass.py 416 
python ../esass.py 416 
