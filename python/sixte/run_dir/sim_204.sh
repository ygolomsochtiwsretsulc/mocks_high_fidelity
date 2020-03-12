#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 204 
python ../pre-process-esass.py 204 
python ../esass.py 204 
