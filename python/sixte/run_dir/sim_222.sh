#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 222 
python ../pre-process-esass.py 222 
python ../esass.py 222 
