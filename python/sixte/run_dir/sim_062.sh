#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 062 
python ../pre-process-esass.py 062 
python ../esass.py 062 
