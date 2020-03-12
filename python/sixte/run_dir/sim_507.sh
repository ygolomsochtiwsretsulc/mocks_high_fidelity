#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 507 
python ../pre-process-esass.py 507 
python ../esass.py 507 
