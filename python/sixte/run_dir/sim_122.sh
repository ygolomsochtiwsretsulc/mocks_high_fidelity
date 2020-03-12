#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 122 
python ../pre-process-esass.py 122 
python ../esass.py 122 
