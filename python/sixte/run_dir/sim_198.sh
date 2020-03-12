#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 198 
python ../pre-process-esass.py 198 
python ../esass.py 198 
