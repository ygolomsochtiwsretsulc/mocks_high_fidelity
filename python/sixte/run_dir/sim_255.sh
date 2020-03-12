#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 255 
python ../pre-process-esass.py 255 
python ../esass.py 255 
