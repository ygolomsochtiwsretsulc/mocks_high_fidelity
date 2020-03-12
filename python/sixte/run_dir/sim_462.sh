#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 462 
python ../pre-process-esass.py 462 
python ../esass.py 462 
