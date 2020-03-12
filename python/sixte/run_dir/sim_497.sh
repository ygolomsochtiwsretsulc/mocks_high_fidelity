#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 497 
python ../pre-process-esass.py 497 
python ../esass.py 497 
