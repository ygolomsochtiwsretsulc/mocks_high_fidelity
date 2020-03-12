#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 111 
python ../pre-process-esass.py 111 
python ../esass.py 111 
