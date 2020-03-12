#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 230 
python ../pre-process-esass.py 230 
python ../esass.py 230 
