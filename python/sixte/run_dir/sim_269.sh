#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 269 
python ../pre-process-esass.py 269 
python ../esass.py 269 
