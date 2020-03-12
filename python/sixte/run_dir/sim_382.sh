#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 382 
python ../pre-process-esass.py 382 
python ../esass.py 382 
