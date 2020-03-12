#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 538 
python ../pre-process-esass.py 538 
python ../esass.py 538 
