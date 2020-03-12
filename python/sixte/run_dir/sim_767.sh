#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 767 
python ../pre-process-esass.py 767 
python ../esass.py 767 
