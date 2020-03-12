#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 712 
python ../pre-process-esass.py 712 
python ../esass.py 712 
