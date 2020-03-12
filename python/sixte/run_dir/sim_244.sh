#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 244 
python ../pre-process-esass.py 244 
python ../esass.py 244 
