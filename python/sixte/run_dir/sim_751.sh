#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 751 
python ../pre-process-esass.py 751 
python ../esass.py 751 
