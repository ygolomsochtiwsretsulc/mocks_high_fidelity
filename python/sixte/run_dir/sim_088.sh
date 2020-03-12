#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 088 
python ../pre-process-esass.py 088 
python ../esass.py 088 
