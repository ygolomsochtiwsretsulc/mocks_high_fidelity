#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 160 
python ../pre-process-esass.py 160 
python ../esass.py 160 
