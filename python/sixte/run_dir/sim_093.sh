#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 093 
python ../pre-process-esass.py 093 
python ../esass.py 093 
