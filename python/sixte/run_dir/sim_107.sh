#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 107 
python ../pre-process-esass.py 107 
python ../esass.py 107 
