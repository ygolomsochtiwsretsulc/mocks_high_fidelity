#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 612 
python ../pre-process-esass.py 612 
python ../esass.py 612 
