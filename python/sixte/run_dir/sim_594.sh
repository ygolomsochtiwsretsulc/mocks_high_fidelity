#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 594 
python ../pre-process-esass.py 594 
python ../esass.py 594 
