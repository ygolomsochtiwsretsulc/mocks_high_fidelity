#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 293 
python ../pre-process-esass.py 293 
python ../esass.py 293 
