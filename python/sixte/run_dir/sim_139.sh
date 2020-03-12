#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 139 
python ../pre-process-esass.py 139 
python ../esass.py 139 
