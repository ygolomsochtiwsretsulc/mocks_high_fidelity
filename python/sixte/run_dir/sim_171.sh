#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 171 
python ../pre-process-esass.py 171 
python ../esass.py 171 
