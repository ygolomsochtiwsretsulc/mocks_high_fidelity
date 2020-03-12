#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 506 
python ../pre-process-esass.py 506 
python ../esass.py 506 
