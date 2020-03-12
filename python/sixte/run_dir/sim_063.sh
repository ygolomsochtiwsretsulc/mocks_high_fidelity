#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 063 
python ../pre-process-esass.py 063 
python ../esass.py 063 
