#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 270 
python ../pre-process-esass.py 270 
python ../esass.py 270 
