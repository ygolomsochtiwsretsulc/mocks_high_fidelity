#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 142 
python ../pre-process-esass.py 142 
python ../esass.py 142 
