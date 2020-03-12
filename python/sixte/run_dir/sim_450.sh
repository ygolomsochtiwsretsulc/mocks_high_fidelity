#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 450 
python ../pre-process-esass.py 450 
python ../esass.py 450 
