#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 728 
python ../pre-process-esass.py 728 
python ../esass.py 728 
