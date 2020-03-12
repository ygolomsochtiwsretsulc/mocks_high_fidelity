#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 453 
python ../pre-process-esass.py 453 
python ../esass.py 453 
