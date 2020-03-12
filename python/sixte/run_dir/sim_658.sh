#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 658 
python ../pre-process-esass.py 658 
python ../esass.py 658 
