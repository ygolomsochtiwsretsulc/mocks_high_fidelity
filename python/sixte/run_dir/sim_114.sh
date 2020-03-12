#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 114 
python ../pre-process-esass.py 114 
python ../esass.py 114 
