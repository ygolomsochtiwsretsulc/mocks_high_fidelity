#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 186 
python ../pre-process-esass.py 186 
python ../esass.py 186 
