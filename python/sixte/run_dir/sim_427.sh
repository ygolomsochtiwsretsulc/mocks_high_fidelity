#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 427 
python ../pre-process-esass.py 427 
python ../esass.py 427 
