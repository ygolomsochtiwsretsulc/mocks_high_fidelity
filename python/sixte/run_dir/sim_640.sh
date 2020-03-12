#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 640 
python ../pre-process-esass.py 640 
python ../esass.py 640 
