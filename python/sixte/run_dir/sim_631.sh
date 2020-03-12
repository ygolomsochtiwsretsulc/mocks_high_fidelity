#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 631 
python ../pre-process-esass.py 631 
python ../esass.py 631 
