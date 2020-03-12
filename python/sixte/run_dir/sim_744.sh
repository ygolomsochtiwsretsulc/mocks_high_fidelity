#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 744 
python ../pre-process-esass.py 744 
python ../esass.py 744 
