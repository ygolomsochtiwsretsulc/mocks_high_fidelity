#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 079 
python ../pre-process-esass.py 079 
python ../esass.py 079 
