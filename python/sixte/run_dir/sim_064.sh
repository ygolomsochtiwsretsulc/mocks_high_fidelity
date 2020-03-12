#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 064 
python ../pre-process-esass.py 064 
python ../esass.py 064 
