#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 710 
python ../pre-process-esass.py 710 
python ../esass.py 710 
