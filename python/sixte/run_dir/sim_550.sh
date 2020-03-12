#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 550 
python ../pre-process-esass.py 550 
python ../esass.py 550 
