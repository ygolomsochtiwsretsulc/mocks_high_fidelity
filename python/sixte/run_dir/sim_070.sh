#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 070 
python ../pre-process-esass.py 070 
python ../esass.py 070 
