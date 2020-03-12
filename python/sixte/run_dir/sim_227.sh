#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 227 
python ../pre-process-esass.py 227 
python ../esass.py 227 
