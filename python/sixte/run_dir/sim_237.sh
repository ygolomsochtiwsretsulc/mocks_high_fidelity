#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 237 
python ../pre-process-esass.py 237 
python ../esass.py 237 
