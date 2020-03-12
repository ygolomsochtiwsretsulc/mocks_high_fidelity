#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 089 
python ../pre-process-esass.py 089 
python ../esass.py 089 
