#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 502 
python ../pre-process-esass.py 502 
python ../esass.py 502 
