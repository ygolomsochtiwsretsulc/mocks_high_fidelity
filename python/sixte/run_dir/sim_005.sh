#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 005 
python ../pre-process-esass.py 005 
python ../esass.py 005 
