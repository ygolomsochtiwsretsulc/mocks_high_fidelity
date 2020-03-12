#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 126 
python ../pre-process-esass.py 126 
python ../esass.py 126 
