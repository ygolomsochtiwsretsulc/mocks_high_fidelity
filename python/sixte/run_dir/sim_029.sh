#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 029 
python ../pre-process-esass.py 029 
python ../esass.py 029 
