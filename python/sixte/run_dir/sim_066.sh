#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 066 
python ../pre-process-esass.py 066 
python ../esass.py 066 
