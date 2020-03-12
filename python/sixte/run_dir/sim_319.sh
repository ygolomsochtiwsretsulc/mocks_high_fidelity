#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 319 
python ../pre-process-esass.py 319 
python ../esass.py 319 
