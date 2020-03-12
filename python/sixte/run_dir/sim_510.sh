#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 510 
python ../pre-process-esass.py 510 
python ../esass.py 510 
