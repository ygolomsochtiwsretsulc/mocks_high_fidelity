#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 423 
python ../pre-process-esass.py 423 
python ../esass.py 423 
