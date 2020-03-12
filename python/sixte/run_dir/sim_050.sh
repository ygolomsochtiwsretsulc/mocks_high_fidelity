#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 050 
python ../pre-process-esass.py 050 
python ../esass.py 050 
