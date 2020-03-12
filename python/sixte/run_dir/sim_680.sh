#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 680 
python ../pre-process-esass.py 680 
python ../esass.py 680 
