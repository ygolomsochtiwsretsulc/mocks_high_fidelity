#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 620 
python ../pre-process-esass.py 620 
python ../esass.py 620 
