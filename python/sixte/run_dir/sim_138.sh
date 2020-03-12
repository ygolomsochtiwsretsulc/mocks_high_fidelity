#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 138 
python ../pre-process-esass.py 138 
python ../esass.py 138 
