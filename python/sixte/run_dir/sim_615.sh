#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 615 
python ../pre-process-esass.py 615 
python ../esass.py 615 
