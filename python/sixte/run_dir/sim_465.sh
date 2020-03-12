#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 465 
python ../pre-process-esass.py 465 
python ../esass.py 465 
