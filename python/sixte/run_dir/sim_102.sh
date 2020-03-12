#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 102 
python ../pre-process-esass.py 102 
python ../esass.py 102 
