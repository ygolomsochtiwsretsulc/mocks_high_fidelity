#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 438 
python ../pre-process-esass.py 438 
python ../esass.py 438 
