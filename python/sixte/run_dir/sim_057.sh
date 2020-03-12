#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 057 
python ../pre-process-esass.py 057 
python ../esass.py 057 
