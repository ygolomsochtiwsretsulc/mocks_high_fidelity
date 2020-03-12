#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 241 
python ../pre-process-esass.py 241 
python ../esass.py 241 
