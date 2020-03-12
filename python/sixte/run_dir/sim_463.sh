#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 463 
python ../pre-process-esass.py 463 
python ../esass.py 463 
