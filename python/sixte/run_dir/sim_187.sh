#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 187 
python ../pre-process-esass.py 187 
python ../esass.py 187 
