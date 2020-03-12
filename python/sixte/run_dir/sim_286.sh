#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 286 
python ../pre-process-esass.py 286 
python ../esass.py 286 
