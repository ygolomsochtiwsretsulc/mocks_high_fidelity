#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 287 
python ../pre-process-esass.py 287 
python ../esass.py 287 
