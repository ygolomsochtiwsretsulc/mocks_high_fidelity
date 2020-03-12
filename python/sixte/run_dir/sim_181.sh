#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 181 
python ../pre-process-esass.py 181 
python ../esass.py 181 
