#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 630 
python ../pre-process-esass.py 630 
python ../esass.py 630 
