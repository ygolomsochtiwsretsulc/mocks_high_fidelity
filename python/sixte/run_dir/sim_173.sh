#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 173 
python ../pre-process-esass.py 173 
python ../esass.py 173 
