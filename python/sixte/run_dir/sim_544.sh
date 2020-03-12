#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 544 
python ../pre-process-esass.py 544 
python ../esass.py 544 
