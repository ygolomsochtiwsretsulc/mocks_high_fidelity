#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 041 
python ../pre-process-esass.py 041 
python ../esass.py 041 
