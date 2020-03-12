#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 432 
python ../pre-process-esass.py 432 
python ../esass.py 432 
