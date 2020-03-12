#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 116 
python ../pre-process-esass.py 116 
python ../esass.py 116 
