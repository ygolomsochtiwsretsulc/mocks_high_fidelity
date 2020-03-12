#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 109 
python ../pre-process-esass.py 109 
python ../esass.py 109 
