#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 168 
python ../pre-process-esass.py 168 
python ../esass.py 168 
