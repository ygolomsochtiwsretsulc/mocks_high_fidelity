#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 129 
python ../pre-process-esass.py 129 
python ../esass.py 129 
