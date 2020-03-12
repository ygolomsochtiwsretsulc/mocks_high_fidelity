#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 113 
python ../pre-process-esass.py 113 
python ../esass.py 113 
