#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 621 
python ../pre-process-esass.py 621 
python ../esass.py 621 
