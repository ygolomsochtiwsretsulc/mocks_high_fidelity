#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 639 
python ../pre-process-esass.py 639 
python ../esass.py 639 
