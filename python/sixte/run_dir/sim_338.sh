#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 338 
python ../pre-process-esass.py 338 
python ../esass.py 338 
