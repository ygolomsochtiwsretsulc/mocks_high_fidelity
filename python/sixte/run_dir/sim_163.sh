#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 163 
python ../pre-process-esass.py 163 
python ../esass.py 163 
