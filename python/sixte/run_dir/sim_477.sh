#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 477 
python ../pre-process-esass.py 477 
python ../esass.py 477 
