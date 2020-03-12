#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 560 
python ../pre-process-esass.py 560 
python ../esass.py 560 
