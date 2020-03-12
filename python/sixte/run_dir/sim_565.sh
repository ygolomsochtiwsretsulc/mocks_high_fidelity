#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 565 
python ../pre-process-esass.py 565 
python ../esass.py 565 
