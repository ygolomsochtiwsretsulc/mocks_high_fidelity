#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 169 
python ../pre-process-esass.py 169 
python ../esass.py 169 
