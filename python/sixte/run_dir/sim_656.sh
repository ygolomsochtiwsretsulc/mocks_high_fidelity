#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 656 
python ../pre-process-esass.py 656 
python ../esass.py 656 
