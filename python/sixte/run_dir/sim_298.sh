#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 298 
python ../pre-process-esass.py 298 
python ../esass.py 298 
