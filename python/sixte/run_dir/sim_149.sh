#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 149 
python ../pre-process-esass.py 149 
python ../esass.py 149 
