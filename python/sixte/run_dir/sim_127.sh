#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 127 
python ../pre-process-esass.py 127 
python ../esass.py 127 
