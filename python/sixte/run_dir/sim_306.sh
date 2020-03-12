#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 306 
python ../pre-process-esass.py 306 
python ../esass.py 306 
