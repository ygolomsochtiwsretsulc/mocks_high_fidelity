#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 702 
python ../pre-process-esass.py 702 
python ../esass.py 702 
