#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 183 
python ../pre-process-esass.py 183 
python ../esass.py 183 
