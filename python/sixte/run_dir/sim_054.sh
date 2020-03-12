#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 054 
python ../pre-process-esass.py 054 
python ../esass.py 054 
