#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 053 
python ../pre-process-esass.py 053 
python ../esass.py 053 
