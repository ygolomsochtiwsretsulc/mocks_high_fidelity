#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 224 
python ../pre-process-esass.py 224 
python ../esass.py 224 
