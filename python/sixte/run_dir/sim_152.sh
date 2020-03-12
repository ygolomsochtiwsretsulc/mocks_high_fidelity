#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 152 
python ../pre-process-esass.py 152 
python ../esass.py 152 
