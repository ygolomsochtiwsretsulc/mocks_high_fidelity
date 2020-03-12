#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 146 
python ../pre-process-esass.py 146 
python ../esass.py 146 
