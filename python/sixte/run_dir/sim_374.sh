#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 374 
python ../pre-process-esass.py 374 
python ../esass.py 374 
