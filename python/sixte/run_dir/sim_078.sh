#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 078 
python ../pre-process-esass.py 078 
python ../esass.py 078 
