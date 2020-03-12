#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 727 
python ../pre-process-esass.py 727 
python ../esass.py 727 
