#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 334 
python ../pre-process-esass.py 334 
python ../esass.py 334 
