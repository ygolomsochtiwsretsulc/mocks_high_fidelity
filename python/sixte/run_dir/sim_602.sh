#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 602 
python ../pre-process-esass.py 602 
python ../esass.py 602 
