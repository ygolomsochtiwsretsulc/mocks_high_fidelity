#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 470 
python ../pre-process-esass.py 470 
python ../esass.py 470 
