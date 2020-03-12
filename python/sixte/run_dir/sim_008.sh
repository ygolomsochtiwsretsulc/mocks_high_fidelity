#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 008 
python ../pre-process-esass.py 008 
python ../esass.py 008 
