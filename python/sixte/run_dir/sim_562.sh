#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 562 
python ../pre-process-esass.py 562 
python ../esass.py 562 
