#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 764 
python ../pre-process-esass.py 764 
python ../esass.py 764 
