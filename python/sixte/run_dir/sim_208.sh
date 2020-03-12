#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 208 
python ../pre-process-esass.py 208 
python ../esass.py 208 
