#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 457 
python ../pre-process-esass.py 457 
python ../esass.py 457 
