#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 275 
python ../pre-process-esass.py 275 
python ../esass.py 275 
